const std = @import("std");
const paper = @import("paper.zig");
pub const Error = paper.Error;
pub const Family = enum(u8) { factual, contradiction, theorem2, theorem3, theorem4, theorem5, repeated_outcome };
pub const Result = struct { bounds: paper.Bound, family: Family, after_consistency: bool = false };
pub const max_cells: usize = 1_048_576;

pub fn checkSize(m: usize, n: usize) Error!void {
    if (m < 2 or n < 2) return error.InvalidDimensions;
    const cells = std.math.mul(usize, m, n) catch return error.WorkspaceLimit;
    if (cells > max_cells) return error.WorkspaceLimit;
}

/// Same formulas and family labels as the preserved R comparator.
/// Column-major matrices; one-based indices; zero denotes no factual condition.
/// Scratch is caller-owned. Only differences is sorted, never an input vector.
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
    selected: []bool,
    differences: []f64,
) Error!Result {
    try checkSize(m, n);
    if (selected.len != m or differences.len < x.len) return error.InvalidDimensions;
    try paper.validateInput(m, n, o, a, x, y, ox, oy, tolerance, px, py);
    var skip: ?usize = null;
    var factual_y = oy;
    for (x, y, 0..) |xi, yi, t| {
        if (xi == ox) {
            if (oy != 0 and yi != oy)
                return .{ .bounds = .{ .lower = 0, .upper = 0 }, .family = .contradiction };
            skip = t;
            factual_y = yi;
        }
    }
    const k = x.len - @intFromBool(skip != null);
    if (k == 0) {
        const mass = if (ox != 0 and factual_y != 0) o[@as(usize, @intCast(ox - 1)) + @as(usize, @intCast(factual_y - 1)) * m] else if (ox != 0) px[@intCast(ox - 1)] else if (factual_y != 0) py[@intCast(factual_y - 1)] else 1;
        return .{ .bounds = .{ .lower = @min(1, @max(0, mass)), .upper = @min(1, @max(0, mass)) }, .family = .factual, .after_consistency = skip != null };
    }
    @memset(selected, false);
    var sum_a: f64 = 0;
    var sum_own: f64 = 0;
    var sum_b: f64 = 0;
    var upper: f64 = 1;
    var subset_count: usize = 0;
    var matching_count: usize = 0;
    for (x, y, 0..) |xi, yi, t| {
        if (skip == t) continue;
        const j: usize = @intCast(xi - 1);
        const cell = j + @as(usize, @intCast(yi - 1)) * m;
        const d = a[cell] - o[cell];
        selected[j] = true;
        sum_a += a[cell];
        sum_own += o[cell];
        sum_b += d + px[j];
        if (ox != 0) {
            upper = @min(upper, d);
        } else if (factual_y == 0 or yi == factual_y) {
            upper = @min(upper, a[cell]);
            differences[subset_count] = d;
            subset_count += 1;
        } else {
            upper = @min(upper, d);
        }
        if (yi == factual_y) matching_count += 1;
    }
    const count: f64 = @floatFromInt(k);
    var lower: f64 = 0;
    var family: Family = undefined;
    if (ox != 0) {
        const mass = if (factual_y == 0) px[@intCast(ox - 1)] else o[@as(usize, @intCast(ox - 1)) + @as(usize, @intCast(factual_y - 1)) * m];
        lower = @max(0, sum_b + mass - count);
        upper = @min(upper, mass);
        family = if (factual_y == 0) .theorem3 else .theorem5;
    } else {
        var mass: f64 = if (factual_y == 0) sum_own else 0;
        for (selected, 0..) |included, j| {
            if (!included) mass += if (factual_y == 0) px[j] else o[j + @as(usize, @intCast(factual_y - 1)) * m];
        }
        if (factual_y != 0) {
            for (x, y) |xi, yi| {
                if (yi == factual_y) mass += o[@as(usize, @intCast(xi - 1)) + @as(usize, @intCast(yi - 1)) * m];
            }
        }
        lower = @max(0, if (factual_y == 0) sum_a - count + 1 else sum_b + mass - count);
        upper = @min(upper, mass);
        for (x, y) |xi, yi| {
            if (factual_y == 0 or yi == factual_y) {
                const j: usize = @intCast(xi - 1);
                const cell = j + @as(usize, @intCast(yi - 1)) * m;
                const b = a[cell] - o[cell] + px[j];
                lower = @max(lower, sum_b - b + o[cell] - count + 1);
            }
        }
        std.mem.sort(f64, differences[0..subset_count], {}, std.sort.asc(f64));
        var prefix: f64 = 0;
        for (differences[0..subset_count], 0..) |d, i| {
            prefix += d;
            if (i > 0) upper = @min(upper, prefix / @as(f64, @floatFromInt(i)));
        }
        family = if (factual_y == 0) .theorem2 else if (matching_count > 1) .repeated_outcome else .theorem4;
    }
    if (!std.math.isFinite(lower) or !std.math.isFinite(upper) or lower > upper + tolerance or lower < -tolerance or upper > 1 + tolerance)
        return error.InconsistentBounds;
    return .{ .bounds = .{ .lower = @min(1, @max(0, lower)), .upper = @min(1, @max(0, upper)) }, .family = family, .after_consistency = skip != null };
}

test "closed form retains distinct-subset upper candidates and theorem labels" {
    const o = [_]f64{ 238.0 / 900.0, 10.0 / 900.0, 147.0 / 900.0, 20.0 / 900.0, 77.0 / 900.0, 72.0 / 900.0, 7.0 / 900.0, 259.0 / 900.0, 70.0 / 900.0 };
    const a = [_]f64{ 80.0 / 300.0, 184.0 / 300.0, 87.0 / 300.0, 7.0 / 300.0, 29.0 / 300.0, 189.0 / 300.0, 213.0 / 300.0, 87.0 / 300.0, 24.0 / 300.0 };
    var px: [3]f64 = undefined;
    var py: [3]f64 = undefined;
    var selected: [3]bool = undefined;
    var differences: [3]f64 = undefined;
    const result = try evaluate(3, 3, &o, &a, &.{ 1, 2, 3 }, &.{ 3, 1, 2 }, 0, 0, 1e-9, &px, &py, &selected, &differences);
    try std.testing.expectApproxEqAbs(@as(f64, 57.0 / 900.0), result.bounds.lower, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 89.0 / 900.0), result.bounds.upper, 1e-12);
    try std.testing.expectEqual(Family.theorem2, result.family);
    try std.testing.expect(!result.after_consistency);
}

test "repeated matching outcomes retain the project-extension label" {
    const o = [_]f64{1.0 / 9.0} ** 9;
    const a = [_]f64{1.0 / 3.0} ** 9;
    var px: [3]f64 = undefined;
    var py: [3]f64 = undefined;
    var selected: [3]bool = undefined;
    var differences: [2]f64 = undefined;
    const result = try evaluate(3, 3, &o, &a, &.{ 1, 2 }, &.{ 1, 1 }, 0, 1, 1e-9, &px, &py, &selected, &differences);
    try std.testing.expectEqual(Family.repeated_outcome, result.family);
    try std.testing.expectApproxEqAbs(@as(f64, 0), result.bounds.lower, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0 / 3.0), result.bounds.upper, 1e-12);
}

test "consistency reduction and contradiction remain distinguishable" {
    const o = [_]f64{ 0.1, 0.2, 0.3, 0.4 };
    const a = [_]f64{ 0.5, 0.3, 0.5, 0.7 };
    var px: [2]f64 = undefined;
    var py: [2]f64 = undefined;
    var selected: [2]bool = undefined;
    var differences: [1]f64 = undefined;
    const result = try evaluate(2, 2, &o, &a, &.{1}, &.{2}, 1, 0, 1e-9, &px, &py, &selected, &differences);
    try std.testing.expectEqual(Family.factual, result.family);
    try std.testing.expect(result.after_consistency);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), result.bounds.lower, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), result.bounds.upper, 1e-12);
    const contradiction = try evaluate(2, 2, &o, &a, &.{1}, &.{2}, 1, 1, 1e-9, &px, &py, &selected, &differences);
    try std.testing.expectEqual(Family.contradiction, contradiction.family);
    try std.testing.expectEqual(@as(f64, 0), contradiction.bounds.upper);
    try std.testing.expect(!contradiction.after_consistency);
}
