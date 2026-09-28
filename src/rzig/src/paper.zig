const std = @import("std");
const scalar = @import("bounds.zig");
pub const Bound = scalar.Bound;
pub const Error = error{ InvalidDimensions, InvalidQuery, InvalidProbability, InvalidTolerance, InvalidNormalization, InconsistentMargins, InconsistentBounds, WorkspaceLimit, Interrupted };
pub const Poll = *const fn () error{Interrupted}!void;
pub const max_cache_states: usize = 1_048_576;

/// The memo is caller-owned and bounded before any allocation or bit shift.
pub fn workspaceSize(m: usize, n: usize, k: usize) Error!usize {
    if (m < 2 or n < 2) return error.InvalidDimensions;
    if (k >= @bitSizeOf(usize) or m >= max_cache_states or n >= max_cache_states)
        return error.WorkspaceLimit;
    const states = std.math.mul(usize, m + 1, n + 1) catch return error.WorkspaceLimit;
    const subsets = @as(usize, 1) << @intCast(k);
    if (states > max_cache_states / subsets) return error.WorkspaceLimit;
    return states * subsets;
}

/// Printed Li-Pearl Theorems 4-11, matching the preserved R recursion.
/// Matrices are column-major. Query labels are one-based; zero is no factual
/// condition. Allocate px, py and memo once before calling this function.
pub fn evaluate(
    m: usize,
    n: usize,
    o: []const f64,
    a: []const f64,
    x: []const i32,
    y: []const i32,
    ox: i32,
    oy: i32,
    tolerance: f64,
    px: []f64,
    py: []f64,
    memo: []Bound,
) Error!Bound {
    return evaluateInterruptible(m, n, o, a, x, y, ox, oy, tolerance, px, py, memo, null);
}

pub fn evaluateInterruptible(
    m: usize,
    n: usize,
    o: []const f64,
    a: []const f64,
    x: []const i32,
    y: []const i32,
    ox: i32,
    oy: i32,
    tolerance: f64,
    px: []f64,
    py: []f64,
    memo: []Bound,
    poll: ?Poll,
) Error!Bound {
    const size = try workspaceSize(m, n, x.len);
    if (memo.len < size) return error.InvalidDimensions;
    try validateInput(m, n, o, a, x, y, ox, oy, tolerance, px, py);
    @memset(memo[0..size], .{ .lower = std.math.nan(f64), .upper = std.math.nan(f64) });
    var engine = Engine{ .m = m, .n = n, .o = o, .a = a, .x = x, .y = y, .px = px, .py = py, .memo = memo[0..size], .tolerance = tolerance, .poll = poll };
    return engine.visit((@as(usize, 1) << @intCast(x.len)) - 1, @intCast(ox), @intCast(oy));
}

/// Shared complete-margin validation, with linear-time duplicate detection.
/// This helper imposes no recursion/memo budget on closed-form consumers.
pub fn validateInput(
    m: usize,
    n: usize,
    o: []const f64,
    a: []const f64,
    x: []const i32,
    y: []const i32,
    ox: i32,
    oy: i32,
    tolerance: f64,
    px: []f64,
    py: []f64,
) Error!void {
    if (m < 2 or n < 2) return error.InvalidDimensions;
    const cells = std.math.mul(usize, m, n) catch return error.InvalidDimensions;
    if (o.len != cells or a.len != cells or px.len != m or py.len != n)
        return error.InvalidDimensions;
    if (!std.math.isFinite(tolerance) or tolerance <= 0) return error.InvalidTolerance;
    if (x.len != y.len or x.len > m or ox < 0 or oy < 0 or ox > m or oy > n)
        return error.InvalidQuery;
    @memset(px, 0);
    for (x, y) |xi, yi| {
        if (xi < 1 or xi > m or yi < 1 or yi > n) return error.InvalidQuery;
        const index: usize = @intCast(xi - 1);
        if (px[index] == 1) return error.InvalidQuery;
        px[index] = 1;
    }
    @memset(px, 0);
    @memset(py, 0);
    var total: f64 = 0;
    for (o, a, 0..) |obs, inter, i| {
        if (!std.math.isFinite(obs) or !std.math.isFinite(inter) or
            obs < -tolerance or inter < -tolerance or obs > 1 + tolerance or inter > 1 + tolerance)
            return error.InvalidProbability;
        px[i % m] += obs;
        py[i / m] += obs;
        total += obs;
    }
    if (@abs(total - 1) > tolerance) return error.InvalidNormalization;
    for (0..m) |j| {
        var sum: f64 = 0;
        for (0..n) |i| {
            const obs = o[j + i * m];
            const inter = a[j + i * m];
            sum += inter;
            if (obs - inter > tolerance or inter - obs - (1 - px[j]) > tolerance)
                return error.InconsistentMargins;
        }
        if (@abs(sum - 1) > tolerance) return error.InvalidNormalization;
    }
}

const Engine = struct {
    m: usize,
    n: usize,
    o: []const f64,
    a: []const f64,
    x: []const i32,
    y: []const i32,
    px: []const f64,
    py: []const f64,
    memo: []Bound,
    tolerance: f64,
    poll: ?Poll,
    visits: usize = 0,

    fn bit(t: usize) usize {
        return @as(usize, 1) << @intCast(t);
    }

    fn oi(self: *const Engine, x: usize, y: usize) f64 {
        return self.o[x - 1 + (y - 1) * self.m];
    }

    fn ai(self: *const Engine, t: usize) f64 {
        return self.a[@as(usize, @intCast(self.x[t] - 1)) + @as(usize, @intCast(self.y[t] - 1)) * self.m];
    }

    fn term(self: *const Engine, mask: usize, x: usize) ?usize {
        for (self.x, 0..) |xi, t| {
            if (mask & bit(t) != 0 and xi == x) return t;
        }
        return null;
    }

    fn single(self: *const Engine, t: usize, ox: usize, oy: usize) Bound {
        const j: usize = @intCast(self.x[t]);
        const i: usize = @intCast(self.y[t]);
        const aji = self.ai(t);
        const oji = self.oi(j, i);
        if (ox == 0 and oy == 0) return .{ .lower = aji, .upper = aji };
        const kind: u8 = if (ox == 0) (if (oy == i) 4 else 5) else (if (oy == 0) 6 else 7);
        if (kind == 5) {
            var partition: f64 = 0;
            for (1..self.m + 1) |p| {
                if (p != j) partition += @max(0, aji + self.oi(p, oy) - 1 + self.px[j - 1] - oji);
            }
            return .{ .lower = @max(0, @max(aji + self.py[oy - 1] - 1, partition)), .upper = @min(aji - oji, self.py[oy - 1] - self.oi(j, oy)) };
        }
        return scalar.single(kind, aji, oji, self.py[i - 1], if (oy == 0) 0 else self.py[oy - 1], self.px[j - 1], if (ox == 0) 0 else self.px[ox - 1], if (ox == 0 or oy == 0) 0 else self.oi(ox, oy), if (oy == 0) 0 else self.oi(j, oy), &.{});
    }

    fn visit(self: *Engine, mask: usize, ox: usize, oy: usize) Error!Bound {
        if (self.visits % 1024 == 0) {
            if (self.poll) |callback| try callback();
        }
        self.visits +%= 1;
        const slot = (mask * (self.m + 1) + ox) * (self.n + 1) + oy;
        if (!std.math.isNan(self.memo[slot].lower)) return self.memo[slot];
        const k = @popCount(mask);
        var result: Bound = undefined;
        if (k == 0) {
            const mass = if (ox != 0 and oy != 0) self.oi(ox, oy) else if (ox != 0) self.px[ox - 1] else if (oy != 0) self.py[oy - 1] else 1;
            result = .{ .lower = mass, .upper = mass };
        } else if (if (ox != 0) self.term(mask, ox) else null) |t| {
            if (oy != 0 and oy != self.y[t]) {
                result = .{ .lower = 0, .upper = 0 };
            } else {
                result = try self.visit(mask ^ bit(t), ox, @intCast(self.y[t]));
            }
        } else if (k == 1) {
            result = self.single(@ctz(mask), ox, oy);
        } else {
            var lower: f64 = 0;
            var upper: f64 = 1;
            var sum_a: f64 = 0;
            for (0..self.x.len) |t| {
                if (mask & bit(t) == 0) continue;
                const margin = self.ai(t);
                sum_a += margin;
                const reduced = try self.visit(mask ^ bit(t), 0, 0);
                upper = @min(upper, @min(margin, reduced.upper));
                if (ox == 0 and oy == 0) {
                    lower = @max(lower, reduced.lower + margin - 1);
                } else {
                    const one = try self.visit(bit(t), ox, oy);
                    lower = @max(lower, reduced.lower + one.lower - 1);
                    upper = @min(upper, one.upper);
                }
            }
            const count: f64 = @floatFromInt(k);
            if (ox == 0 and oy == 0) {
                lower = @max(lower, sum_a - count + 1);
            } else {
                const mass = if (ox != 0 and oy != 0) self.oi(ox, oy) else if (ox != 0) self.px[ox - 1] else self.py[oy - 1];
                lower = @max(lower, mass - count + sum_a);
                upper = @min(upper, mass);
            }
            if (ox == 0) {
                var by_lower: f64 = 0;
                var by_upper: f64 = 0;
                for (1..self.m + 1) |p| {
                    const by = if (self.term(mask, p)) |t| value: {
                        if (oy != 0 and self.y[t] != oy) break :value Bound{ .lower = 0, .upper = 0 };
                        break :value try self.visit(mask ^ bit(t), p, @intCast(self.y[t]));
                    } else try self.visit(mask, p, oy);
                    by_lower += by.lower;
                    by_upper += by.upper;
                }
                lower = @max(lower, by_lower);
                upper = @min(upper, by_upper);
            }
            result = .{ .lower = lower, .upper = upper };
        }
        if (!std.math.isFinite(result.lower) or !std.math.isFinite(result.upper) or
            result.lower > result.upper + self.tolerance or result.lower < -self.tolerance or result.upper > 1 + self.tolerance)
            return error.InconsistentBounds;
        result.lower = @min(1, @max(0, result.lower));
        result.upper = @min(1, @max(0, result.upper));
        self.memo[slot] = result;
        return result;
    }
};

test "recursive conjunction reproduces the published three-treatment example" {
    const o = [_]f64{ 238.0 / 900.0, 10.0 / 900.0, 147.0 / 900.0, 20.0 / 900.0, 77.0 / 900.0, 72.0 / 900.0, 7.0 / 900.0, 259.0 / 900.0, 70.0 / 900.0 };
    const a = [_]f64{ 80.0 / 300.0, 184.0 / 300.0, 87.0 / 300.0, 7.0 / 300.0, 29.0 / 300.0, 189.0 / 300.0, 213.0 / 300.0, 87.0 / 300.0, 24.0 / 300.0 };
    var px: [3]f64 = undefined;
    var py: [3]f64 = undefined;
    var scratch: [128]Bound = undefined;
    const answer = try evaluate(3, 3, &o, &a, &.{ 1, 2, 3 }, &.{ 3, 1, 2 }, 0, 0, 1e-9, &px, &py, &scratch);
    try std.testing.expectApproxEqAbs(@as(f64, 0), answer.lower, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 89.0 / 900.0), answer.upper, 1e-12);
}

test "all recursive factual-condition families preserve the R oracle endpoints" {
    const o = [_]f64{ 238.0 / 900.0, 10.0 / 900.0, 147.0 / 900.0, 20.0 / 900.0, 77.0 / 900.0, 72.0 / 900.0, 7.0 / 900.0, 259.0 / 900.0, 70.0 / 900.0 };
    const a = [_]f64{ 80.0 / 300.0, 184.0 / 300.0, 87.0 / 300.0, 7.0 / 300.0, 29.0 / 300.0, 189.0 / 300.0, 213.0 / 300.0, 87.0 / 300.0, 24.0 / 300.0 };
    var px: [3]f64 = undefined;
    var py: [3]f64 = undefined;
    var scratch: [64]Bound = undefined;
    const conditions = [_][2]i32{ .{ 0, 0 }, .{ 3, 0 }, .{ 0, 1 }, .{ 3, 1 }, .{ 1, 1 }, .{ 1, 3 } };
    const expected = [_][2]f64{ .{ 291.0 / 900.0, 306.0 / 900.0 }, .{ 16.0 / 900.0, 289.0 / 900.0 }, .{ 7.0 / 900.0, 157.0 / 900.0 }, .{ 0, 147.0 / 900.0 }, .{ 0, 0 }, .{ 0, 7.0 / 900.0 } };
    for (conditions, expected) |condition, want| {
        const answer = try evaluate(3, 3, &o, &a, &.{ 1, 2 }, &.{ 3, 1 }, condition[0], condition[1], 1e-9, &px, &py, &scratch);
        try std.testing.expectApproxEqAbs(want[0], answer.lower, 1e-12);
        try std.testing.expectApproxEqAbs(want[1], answer.upper, 1e-12);
    }
}

test "workspace limits reject overflow before allocation" {
    try std.testing.expectEqual(@as(usize, 128), try workspaceSize(3, 3, 3));
    try std.testing.expectError(error.WorkspaceLimit, workspaceSize(100, 100, 20));
    try std.testing.expectError(error.WorkspaceLimit, workspaceSize(3, 3, 64));
    try std.testing.expectError(error.InvalidDimensions, workspaceSize(1, 3, 1));
}

test "interruption is propagated without a partial result" {
    const interrupt = struct {
        fn poll() error{Interrupted}!void {
            return error.Interrupted;
        }
    };
    var px: [2]f64 = undefined;
    var py: [2]f64 = undefined;
    var scratch: [18]Bound = undefined;
    try std.testing.expectError(error.Interrupted, evaluateInterruptible(2, 2, &.{ 0.1, 0.2, 0.3, 0.4 }, &.{ 0.5, 0.3, 0.5, 0.7 }, &.{1}, &.{2}, 0, 0, 1e-9, &px, &py, &scratch, interrupt.poll));
}
