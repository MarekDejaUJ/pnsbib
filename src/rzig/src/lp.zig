const std = @import("std");

pub const Error = error{ InvalidDimensions, InvalidInput, InvalidOptions, WorkspaceLimit, IterationLimit, NumericalFailure, Interrupted };
pub const Poll = *const fn () error{Interrupted}!void;
pub const Status = enum { optimal, infeasible, unbounded };
pub const Options = struct {
    pivot_tolerance: f64 = 1e-10,
    feasibility_tolerance: f64 = 1e-8,
    max_iterations: usize = 200_000,
};
pub const Result = struct {
    status: Status,
    objective: f64 = std.math.nan(f64),
    iterations: usize = 0,
    max_primal_violation: f64 = 0,
    max_dual_violation: f64 = 0,
    duality_gap: f64 = 0,
};
pub const max_tableau_cells: usize = 16_777_216;

pub fn workspaceSize(m: usize, n: usize) Error!usize {
    if (n == 0) return error.InvalidDimensions;
    const aux = std.math.mul(usize, m, 2) catch return error.WorkspaceLimit;
    const columns = std.math.add(usize, n, aux) catch return error.WorkspaceLimit;
    const width = std.math.add(usize, columns, 1) catch return error.WorkspaceLimit;
    const rows = std.math.add(usize, m, 1) catch return error.WorkspaceLimit;
    const size = std.math.mul(usize, rows, width) catch return error.WorkspaceLimit;
    const lex_size = std.math.mul(usize, m, m) catch return error.WorkspaceLimit;
    const total = std.math.add(usize, size, lex_size) catch return error.WorkspaceLimit;
    if (total > max_tableau_cells) return error.WorkspaceLimit;
    return size;
}

/// All memory is provided at setup. Pivot iterations do not allocate.
pub const Workspace = struct {
    tableau: []f64,
    basis: []usize,
    basic: []bool,
    active: []bool,
    scale: []f64,
    lex: []f64,

    pub fn initAlloc(alloc: std.mem.Allocator, m: usize, n: usize) !Workspace {
        const size = try workspaceSize(m, n);
        const tableau = try alloc.alloc(f64, size);
        errdefer alloc.free(tableau);
        const basis = try alloc.alloc(usize, m);
        errdefer alloc.free(basis);
        const basic = try alloc.alloc(bool, n + 2 * m);
        errdefer alloc.free(basic);
        const active = try alloc.alloc(bool, m);
        errdefer alloc.free(active);
        const scale = try alloc.alloc(f64, m);
        errdefer alloc.free(scale);
        const lex = try alloc.alloc(f64, m * m);
        return .{ .tableau = tableau, .basis = basis, .basic = basic, .active = active, .scale = scale, .lex = lex };
    }

    pub fn deinit(self: *Workspace, alloc: std.mem.Allocator) void {
        alloc.free(self.tableau);
        alloc.free(self.basis);
        alloc.free(self.basic);
        alloc.free(self.active);
        alloc.free(self.scale);
        alloc.free(self.lex);
    }
};

/// Continuous LP with nonnegative variables. A is column-major; directions
/// -1, 0, +1 mean <=, =, >=. For optimal results dual has the original min/max
/// convention; for infeasibility it contains a separating Farkas vector.
/// An unbounded result supplies a feasible primal point and an improving ray.
pub fn solve(
    m: usize,
    n: usize,
    a: []const f64,
    b: []const f64,
    directions: []const i32,
    c: []const f64,
    maximize: bool,
    options: Options,
    primal: []f64,
    dual: []f64,
    ray: []f64,
    scratch: *Workspace,
    poll: ?Poll,
) Error!Result {
    const size = try workspaceSize(m, n);
    if (a.len != m * n or b.len != m or directions.len != m or c.len != n or
        primal.len != n or dual.len != m or ray.len != n or scratch.tableau.len != size or
        scratch.basis.len != m or scratch.basic.len != n + 2 * m or scratch.active.len != m or scratch.scale.len != m or scratch.lex.len != m * m)
        return error.InvalidDimensions;
    if (!std.math.isFinite(options.pivot_tolerance) or !std.math.isFinite(options.feasibility_tolerance) or
        options.pivot_tolerance <= 0 or options.feasibility_tolerance < options.pivot_tolerance or options.max_iterations == 0)
        return error.InvalidOptions;
    for (a) |v| if (!std.math.isFinite(v)) return error.InvalidInput;
    for (b) |v| if (!std.math.isFinite(v)) return error.InvalidInput;
    for (c) |v| if (!std.math.isFinite(v)) return error.InvalidInput;
    for (directions) |d| if (d < -1 or d > 1) return error.InvalidInput;
    if (poll) |callback| try callback();
    @memset(primal, 0);
    @memset(dual, 0);
    @memset(ray, 0);
    @memset(scratch.tableau, 0);
    @memset(scratch.basic, false);
    @memset(scratch.active, true);
    var engine = Engine{ .m = m, .n = n, .width = n + 2 * m + 1, .w = scratch, .options = options, .poll = poll };
    const rhs = engine.width - 1;
    for (0..m) |i| {
        var scale: f64 = 0;
        for (0..n) |j| scale = @max(scale, @abs(a[i + j * m]));
        if (scale == 0) scale = 1;
        const sign: f64 = if (b[i] < 0) -1 else 1;
        scratch.scale[i] = sign / scale;
        if (!std.math.isFinite(scratch.scale[i])) return error.NumericalFailure;
        for (0..n) |j| engine.cell(i, j).* = a[i + j * m] / scale * sign;
        engine.cell(i, rhs).* = @abs(b[i]) / scale;
        if (!std.math.isFinite(engine.cell(i, rhs).*)) return error.NumericalFailure;
        const d = directions[i] * @as(i32, if (sign < 0) -1 else 1);
        if (d != 0) engine.cell(i, n + i).* = if (d == -1) 1 else -1;
        // Retain all artificial identity columns for original-row dual recovery,
        // including rows initially using a slack basis.
        engine.cell(i, n + m + i).* = 1;
        scratch.basis[i] = if (d == -1) n + i else n + m + i;
        scratch.basic[scratch.basis[i]] = true;
        engine.cell(m, n + m + i).* = 1;
    }
    try engine.canonicalize();
    if (try engine.simplex(rhs, true) != null) return error.NumericalFailure;
    const phase_one = engine.cell(m, rhs).*;
    if (phase_one < -options.feasibility_tolerance) {
        for (0..m) |i| dual[i] = (engine.cell(m, n + m + i).* - 1) * scratch.scale[i];
        var result = Result{ .status = .infeasible, .iterations = engine.iterations };
        try checkDual(m, n, a, b, directions, c, 0, dual, options.feasibility_tolerance, &result);
        var certificate: f64 = 0;
        var magnitude: f64 = 0;
        for (b, dual) |bi, yi| {
            certificate += bi * yi;
            magnitude += @abs(bi * yi);
        }
        if (!std.math.isFinite(certificate) or !std.math.isFinite(magnitude) or
            certificate >= -options.feasibility_tolerance * (1 + magnitude)) return error.NumericalFailure;
        return result;
    }
    // A zero artificial basic variable is removed by a non-artificial pivot;
    // a row with no such coefficient is redundant in the original variables.
    for (0..m) |i| {
        if (scratch.basis[i] < n + m) continue;
        if (@abs(engine.cell(i, rhs).*) > options.feasibility_tolerance) return error.NumericalFailure;
        var entering: ?usize = null;
        for (0..n + m) |j| {
            if (!scratch.basic[j] and @abs(engine.cell(i, j).*) > options.pivot_tolerance) {
                entering = j;
                break;
            }
        }
        if (entering) |j| {
            try engine.pivot(i, j);
        } else {
            scratch.active[i] = false;
        }
    }
    @memset(scratch.tableau[m * engine.width ..], 0);
    const sign: f64 = if (maximize) 1 else -1;
    for (c, 0..) |v, j| engine.cell(m, j).* = -sign * v;
    try engine.canonicalize();
    const unbounded_column = try engine.simplex(n + m, false);
    for (0..m) |i| {
        if (scratch.active[i] and scratch.basis[i] < n)
            primal[scratch.basis[i]] = engine.cell(i, rhs).*;
    }
    var result = Result{ .status = if (unbounded_column == null) .optimal else .unbounded, .iterations = engine.iterations };
    try checkPrimal(m, n, a, b, directions, primal, false, options.feasibility_tolerance, &result);
    if (unbounded_column) |entering| {
        if (entering < n) ray[entering] = 1;
        for (0..m) |i| {
            if (scratch.active[i] and scratch.basis[i] < n)
                ray[scratch.basis[i]] = -engine.cell(i, entering).*;
        }
        try checkPrimal(m, n, a, b, directions, ray, true, options.feasibility_tolerance, &result);
        var improvement: f64 = 0;
        for (c, ray) |ci, ri| improvement += sign * ci * ri;
        if (!std.math.isFinite(improvement) or improvement <= options.pivot_tolerance) return error.NumericalFailure;
        return result;
    }
    for (0..m) |i| dual[i] = engine.cell(m, n + m + i).* * scratch.scale[i];
    try checkDual(m, n, a, b, directions, c, sign, dual, options.feasibility_tolerance, &result);
    var primal_value: f64 = 0;
    var dual_value: f64 = 0;
    for (c, primal) |ci, xi| primal_value += sign * ci * xi;
    for (b, dual) |bi, yi| dual_value += bi * yi;
    result.duality_gap = @abs(primal_value - dual_value);
    if (!std.math.isFinite(primal_value) or !std.math.isFinite(dual_value) or
        result.duality_gap > options.feasibility_tolerance * (1 + @abs(primal_value) + @abs(dual_value)))
        return error.NumericalFailure;
    result.objective = sign * primal_value;
    for (dual) |*v| v.* *= sign;
    if (poll) |callback| try callback();
    return result;
}

const Engine = struct {
    m: usize,
    n: usize,
    width: usize,
    w: *Workspace,
    options: Options,
    poll: ?Poll,
    iterations: usize = 0,

    fn cell(self: *Engine, i: usize, j: usize) *f64 {
        return &self.w.tableau[i * self.width + j];
    }

    fn canonicalize(self: *Engine) Error!void {
        for (0..self.m) |i| {
            if (!self.w.active[i]) continue;
            const coefficient = self.cell(self.m, self.w.basis[i]).*;
            if (coefficient == 0) continue;
            for (0..self.width) |j| {
                const v = self.cell(self.m, j);
                v.* -= coefficient * self.cell(i, j).*;
                if (!std.math.isFinite(v.*)) return error.NumericalFailure;
            }
            self.cell(self.m, self.w.basis[i]).* = 0;
        }
    }

    fn pivot(self: *Engine, leaving: usize, entering: usize) Error!void {
        if (self.iterations >= self.options.max_iterations) return error.IterationLimit;
        if (self.iterations % 32 == 0) {
            if (self.poll) |callback| try callback();
        }
        self.iterations += 1;
        const divisor = self.cell(leaving, entering).*;
        if (!std.math.isFinite(divisor) or @abs(divisor) <= self.options.pivot_tolerance) return error.NumericalFailure;
        for (0..self.width) |j| {
            self.cell(leaving, j).* /= divisor;
            if (!std.math.isFinite(self.cell(leaving, j).*)) return error.NumericalFailure;
        }
        self.cell(leaving, entering).* = 1;
        for (self.w.lex[leaving * self.m ..][0..self.m]) |*v| {
            v.* /= divisor;
            if (!std.math.isFinite(v.*)) return error.NumericalFailure;
        }
        for (0..self.m + 1) |i| {
            if (i == leaving or (i < self.m and !self.w.active[i])) continue;
            const coefficient = self.cell(i, entering).*;
            if (coefficient == 0) continue;
            for (0..self.width) |j| {
                const v = self.cell(i, j);
                v.* -= coefficient * self.cell(leaving, j).*;
                if (!std.math.isFinite(v.*)) return error.NumericalFailure;
            }
            if (i < self.m) for (0..self.m) |j| {
                const v = &self.w.lex[i * self.m + j];
                v.* -= coefficient * self.w.lex[leaving * self.m + j];
                if (!std.math.isFinite(v.*)) return error.NumericalFailure;
            };
            self.cell(i, entering).* = 0;
        }
        self.w.basic[self.w.basis[leaving]] = false;
        self.w.basis[leaving] = entering;
        self.w.basic[entering] = true;
    }

    /// Compare only symbolic perturbations after an equal ordinary ratio.
    /// Separate scratch preserves original identity columns for certificates.
    fn lexLess(self: *Engine, left: usize, right: usize, column: usize) Error!bool {
        const left_pivot = self.cell(left, column).*;
        const right_pivot = self.cell(right, column).*;
        for (0..self.m) |j| {
            const lhs = self.w.lex[left * self.m + j] / left_pivot;
            const rhs = self.w.lex[right * self.m + j] / right_pivot;
            if (!std.math.isFinite(lhs) or !std.math.isFinite(rhs)) return error.NumericalFailure;
            if (lhs != rhs) return lhs < rhs;
        }
        return error.NumericalFailure;
    }

    /// Returns the entering column only when no leaving row exists.
    fn simplex(self: *Engine, allowed_columns: usize, feasibility_phase: bool) Error!?usize {
        // Each phase starts from its own feasible canonical system, including
        // phase II after removal of zero artificial basics/redundant rows.
        @memset(self.w.lex, 0);
        for (0..self.m) |i| self.w.lex[i * self.m + i] = 1;
        while (true) {
            // Phase I maximizes minus the nonnegative artificial mass, whose
            // known upper bound is zero. Once that bound is attained within
            // the existing feasibility tolerance, further degenerate pivots
            // only amplify cancellation in redundant/paired constraint rows.
            // Require both actual basic artificial mass and objective agreement;
            // phase II still verifies original-matrix primal/dual certificates.
            if (feasibility_phase) {
                var artificial_mass: f64 = 0;
                for (0..self.m) |i| {
                    if (self.w.active[i] and self.w.basis[i] >= self.n + self.m)
                        artificial_mass += @abs(self.cell(i, self.width - 1).*);
                }
                if (artificial_mass <= self.options.feasibility_tolerance and
                    @abs(self.cell(self.m, self.width - 1).*) <= self.options.feasibility_tolerance)
                    return null;
            }
            var entering: ?usize = null;
            var reduced_cost = -self.options.pivot_tolerance;
            for (0..allowed_columns) |j| {
                if (!self.w.basic[j] and self.cell(self.m, j).* < reduced_cost) {
                    entering = j;
                    reduced_cost = self.cell(self.m, j).*;
                }
            }
            const column = entering orelse return null;
            var leaving: ?usize = null;
            var ratio = std.math.inf(f64);
            for (0..self.m) |i| {
                if (!self.w.active[i]) continue;
                const coefficient = self.cell(i, column).*;
                if (coefficient <= self.options.pivot_tolerance) continue;
                const rhs = self.cell(i, self.width - 1).*;
                if (rhs < -self.options.feasibility_tolerance * (1 + @abs(rhs))) return error.NumericalFailure;
                const candidate = @max(0, rhs) / coefficient;
                if (!std.math.isFinite(candidate)) return error.NumericalFailure;
                if (leaving == null or candidate < ratio or
                    (candidate == ratio and try self.lexLess(i, leaving.?, column)))
                {
                    leaving = i;
                    ratio = candidate;
                }
            }
            if (leaving) |row| {
                // Safeguarded two-pass ratio selection: retain the strict
                // lexicographic choice unless it would divide by a pivot tiny
                // relative to another admissible one. The relaxed step uses
                // the smaller, unchanged pivot tolerance, not a widened
                // feasibility tolerance. Original-matrix certificates remain
                // mandatory before any mathematical status is returned.
                var relaxed = std.math.inf(f64);
                for (0..self.m) |i| {
                    if (!self.w.active[i]) continue;
                    const coefficient = self.cell(i, column).*;
                    if (coefficient <= self.options.pivot_tolerance) continue;
                    const rhs = self.cell(i, self.width - 1).*;
                    relaxed = @min(relaxed, (@max(0, rhs) +
                        self.options.pivot_tolerance * (1 + @abs(rhs))) / coefficient);
                }
                var stable_row = row;
                var largest = self.cell(row, column).*;
                for (0..self.m) |i| {
                    if (!self.w.active[i]) continue;
                    const coefficient = self.cell(i, column).*;
                    if (coefficient <= largest) continue;
                    const step = @max(0, self.cell(i, self.width - 1).*) / coefficient;
                    if (step <= relaxed) {
                        stable_row = i;
                        largest = coefficient;
                    }
                }
                if (self.cell(row, column).* < self.options.pivot_tolerance * largest)
                    leaving = stable_row;
                try self.pivot(leaving.?, column);
            } else return column;
        }
    }
};

fn checkPrimal(m: usize, n: usize, a: []const f64, b: []const f64, directions: []const i32, x: []const f64, recession: bool, tolerance: f64, result: *Result) Error!void {
    for (x) |xi| {
        if (!std.math.isFinite(xi) or xi < -tolerance * (1 + @abs(xi))) return error.NumericalFailure;
        result.max_primal_violation = @max(result.max_primal_violation, -xi);
    }
    for (0..m) |i| {
        const rhs: f64 = if (recession) 0 else b[i];
        var value: f64 = 0;
        var magnitude: f64 = @abs(rhs);
        for (0..n) |j| {
            const term = a[i + j * m] * x[j];
            value += term;
            magnitude += @abs(term);
        }
        const difference = value - rhs;
        const violation = if (directions[i] == 0) @abs(difference) else if (directions[i] == -1) difference else -difference;
        if (!std.math.isFinite(value) or !std.math.isFinite(magnitude) or violation > tolerance * (1 + magnitude)) return error.NumericalFailure;
        result.max_primal_violation = @max(result.max_primal_violation, violation);
    }
}

fn checkDual(m: usize, n: usize, a: []const f64, b: []const f64, directions: []const i32, c: []const f64, sign: f64, dual: []const f64, tolerance: f64, result: *Result) Error!void {
    _ = b;
    for (dual, directions) |yi, direction| {
        const violation = @as(f64, @floatFromInt(direction)) * yi;
        if (!std.math.isFinite(yi) or violation > tolerance * (1 + @abs(yi))) return error.NumericalFailure;
        result.max_dual_violation = @max(result.max_dual_violation, violation);
    }
    for (0..n) |j| {
        var value: f64 = 0;
        var magnitude: f64 = @abs(sign * c[j]);
        for (0..m) |i| {
            const term = a[i + j * m] * dual[i];
            value += term;
            magnitude += @abs(term);
        }
        const violation = sign * c[j] - value;
        if (!std.math.isFinite(value) or !std.math.isFinite(magnitude) or violation > tolerance * (1 + magnitude)) return error.NumericalFailure;
        result.max_dual_violation = @max(result.max_dual_violation, violation);
    }
}

const Fixture = struct {
    result: Result,
    primal: []f64,
    dual: []f64,
    ray: []f64,
    fn deinit(self: *Fixture) void {
        std.testing.allocator.free(self.primal);
        std.testing.allocator.free(self.dual);
        std.testing.allocator.free(self.ray);
    }
};

fn fixture(m: usize, n: usize, a: []const f64, b: []const f64, directions: []const i32, c: []const f64, maximize: bool, options: Options, poll: ?Poll) !Fixture {
    const alloc = std.testing.allocator;
    var workspace = try Workspace.initAlloc(alloc, m, n);
    defer workspace.deinit(alloc);
    const primal = try alloc.alloc(f64, n);
    errdefer alloc.free(primal);
    const dual = try alloc.alloc(f64, m);
    errdefer alloc.free(dual);
    const ray = try alloc.alloc(f64, n);
    errdefer alloc.free(ray);
    const result = try solve(m, n, a, b, directions, c, maximize, options, primal, dual, ray, &workspace, poll);
    return .{ .result = result, .primal = primal, .dual = dual, .ray = ray };
}

test "mixed directions signed right hand sides and redundant equalities" {
    const a = [_]f64{ 1, 1, 1, -1, 1, 0, 0, -1 };
    const b = [_]f64{ 3, 1, 2, -3 };
    const d = [_]i32{ 0, 1, -1, 0 };
    for ([_]bool{ true, false }, [_]f64{ 5, 4 }) |maximize, optimum| {
        var answer = try fixture(4, 2, &a, &b, &d, &.{ 2, 1 }, maximize, .{}, null);
        defer answer.deinit();
        try std.testing.expectEqual(Status.optimal, answer.result.status);
        try std.testing.expectApproxEqAbs(optimum, answer.result.objective, 1e-9);
        try std.testing.expectApproxEqAbs(optimum - 3, answer.primal[0], 1e-9);
    }
}

test "bounded maximum supplies primal and dual witnesses" {
    var answer = try fixture(3, 2, &.{ 1, 1, 0, 1, 0, 1 }, &.{ 4, 2, 3 }, &.{ -1, -1, -1 }, &.{ 3, 2 }, true, .{}, null);
    defer answer.deinit();
    try std.testing.expectEqual(Status.optimal, answer.result.status);
    try std.testing.expectApproxEqAbs(10, answer.result.objective, 1e-9);
    try std.testing.expectApproxEqAbs(2, answer.primal[0], 1e-9);
    try std.testing.expectApproxEqAbs(2, answer.primal[1], 1e-9);
    try std.testing.expectApproxEqAbs(2, answer.dual[0], 1e-9);
    try std.testing.expectApproxEqAbs(1, answer.dual[1], 1e-9);
    try std.testing.expect(answer.result.duality_gap < 1e-9);
}

test "infeasibility has a separating Farkas vector" {
    var answer = try fixture(2, 1, &.{ 1, 1 }, &.{ 0, 1 }, &.{ -1, 1 }, &.{0}, true, .{}, null);
    defer answer.deinit();
    try std.testing.expectEqual(Status.infeasible, answer.result.status);
    try std.testing.expect(answer.dual[0] >= 0 and answer.dual[1] <= 0);
    try std.testing.expect(answer.dual[0] + answer.dual[1] >= -1e-9);
    try std.testing.expect(answer.dual[1] < -1e-8);
}

test "unboundedness has a feasible point and improving ray" {
    var answer = try fixture(1, 2, &.{ 1, -1 }, &.{0}, &.{1}, &.{ 1, 0 }, true, .{}, null);
    defer answer.deinit();
    try std.testing.expectEqual(Status.unbounded, answer.result.status);
    try std.testing.expect(answer.ray[0] > 0 and answer.ray[1] >= 0);
    try std.testing.expect(answer.ray[0] - answer.ray[1] >= -1e-9);
}

test "zero rows and zero objective remain feasible" {
    var answer = try fixture(3, 1, &.{ 0, 0, 0 }, &.{ 0, 0, -1 }, &.{ 0, -1, 1 }, &.{0}, true, .{}, null);
    defer answer.deinit();
    try std.testing.expectEqual(Status.optimal, answer.result.status);
    try std.testing.expectEqual(@as(f64, 0), answer.result.objective);
    var impossible = try fixture(1, 1, &.{0}, &.{1}, &.{0}, &.{0}, true, .{}, null);
    defer impossible.deinit();
    try std.testing.expectEqual(Status.infeasible, impossible.result.status);
}

test "empty constraint set handles both optimal and unbounded objectives" {
    var bounded = try fixture(0, 2, &.{}, &.{}, &.{}, &.{ -1, 0 }, true, .{}, null);
    defer bounded.deinit();
    try std.testing.expectEqual(Status.optimal, bounded.result.status);
    var unbounded = try fixture(0, 2, &.{}, &.{}, &.{}, &.{ 1, 0 }, true, .{}, null);
    defer unbounded.deinit();
    try std.testing.expectEqual(Status.unbounded, unbounded.result.status);
}

test "lexicographic choices terminate the classical cycling fixture" {
    var answer = try fixture(3, 4, &.{ 0.5, 0.5, 1, -5.5, -1.5, 0, -2.5, -0.5, 0, 9, 1, 0 }, &.{ 0, 0, 1 }, &.{ -1, -1, -1 }, &.{ 10, -57, -9, -24 }, true, .{}, null);
    defer answer.deinit();
    try std.testing.expectEqual(Status.optimal, answer.result.status);
    try std.testing.expectApproxEqAbs(1, answer.result.objective, 1e-9);
    try std.testing.expect(answer.result.iterations < 100);
}

test "generic auxiliary variables are not clamped to probabilities" {
    var answer = try fixture(2, 1, &.{ 1e-4, 1e4 }, &.{ 100, 10000010000 }, &.{ 1, -1 }, &.{1}, true, .{}, null);
    defer answer.deinit();
    try std.testing.expectEqual(Status.optimal, answer.result.status);
    try std.testing.expectApproxEqAbs(1000001, answer.primal[0], 1e-6);
}

test "malformed dimensions values directions and tolerances fail explicitly" {
    try std.testing.expectError(error.InvalidDimensions, fixture(1, 2, &.{1}, &.{1}, &.{0}, &.{ 1, 0 }, true, .{}, null));
    try std.testing.expectError(error.InvalidInput, fixture(1, 1, &.{std.math.nan(f64)}, &.{1}, &.{0}, &.{1}, true, .{}, null));
    try std.testing.expectError(error.InvalidInput, fixture(1, 1, &.{1}, &.{1}, &.{2}, &.{1}, true, .{}, null));
    try std.testing.expectError(error.InvalidOptions, fixture(1, 1, &.{1}, &.{1}, &.{0}, &.{1}, true, .{ .pivot_tolerance = 0 }, null));
    try std.testing.expectError(error.WorkspaceLimit, workspaceSize(std.math.maxInt(usize), 2));
    try std.testing.expectError(error.InvalidDimensions, workspaceSize(1, 0));
}

test "iteration limits and interruption are not mathematical infeasibility" {
    try std.testing.expectError(error.IterationLimit, fixture(3, 2, &.{ 1, 1, 0, 1, 0, 1 }, &.{ 4, 2, 3 }, &.{ -1, -1, -1 }, &.{ 3, 2 }, true, .{ .max_iterations = 1 }, null));
    const Stop = struct {
        fn poll() error{Interrupted}!void {
            return error.Interrupted;
        }
    };
    try std.testing.expectError(error.Interrupted, fixture(1, 1, &.{1}, &.{1}, &.{0}, &.{1}, true, .{}, Stop.poll));
}

test "lexicographic ratio ties compare transformed perturbation rows" {
    var workspace = try Workspace.initAlloc(std.testing.allocator, 2, 1);
    defer workspace.deinit(std.testing.allocator);
    @memset(workspace.tableau, 0);
    @memset(workspace.lex, 0);
    var engine = Engine{ .m = 2, .n = 1, .width = 6, .w = &workspace, .options = .{}, .poll = null };
    engine.cell(0, 0).* = 2;
    engine.cell(1, 0).* = 1;
    workspace.lex[0] = 1;
    workspace.lex[3] = 1;
    try std.testing.expect(try engine.lexLess(1, 0, 0));
    try std.testing.expect(!try engine.lexLess(0, 1, 0));
    workspace.lex[2] = 2;
    try std.testing.expect(!try engine.lexLess(1, 0, 0));
}

test "phase one stops at zero artificial mass before degenerate pivots" {
    var w = try Workspace.initAlloc(std.testing.allocator, 1, 1);
    defer w.deinit(std.testing.allocator);
    @memset(w.tableau, 0);
    @memset(w.basic, false);
    @memset(w.active, true);
    w.basis[0] = 2;
    w.basic[2] = true;
    var engine = Engine{ .m = 1, .n = 1, .width = 4, .w = &w, .options = .{}, .poll = null };
    // x - slack + artificial = 0 already attains phase I's upper bound zero.
    // Its negative reduced cost only permits a degenerate zero-length pivot.
    engine.cell(0, 0).* = 1;
    engine.cell(0, 1).* = -1;
    engine.cell(0, 2).* = 1;
    engine.cell(1, 0).* = -1;
    engine.cell(1, 1).* = 1;
    try std.testing.expectEqual(@as(?usize, null), try engine.simplex(3, true));
    try std.testing.expectEqual(@as(usize, 0), engine.iterations);
}

test "ratio safeguard avoids a tiny pivot with a stable admissible alternative" {
    var w = try Workspace.initAlloc(std.testing.allocator, 2, 1);
    defer w.deinit(std.testing.allocator);
    @memset(w.tableau, 0);
    @memset(w.basic, false);
    @memset(w.active, true);
    w.basis[0] = 1;
    w.basis[1] = 2;
    w.basic[1] = true;
    w.basic[2] = true;
    var engine = Engine{ .m = 2, .n = 1, .width = 6, .w = &w, .options = .{}, .poll = null };
    engine.cell(0, 0).* = 1.5e-10;
    engine.cell(0, 1).* = 1;
    engine.cell(1, 0).* = 500;
    engine.cell(1, 2).* = 1;
    engine.cell(1, 5).* = 1e-12;
    engine.cell(2, 0).* = -1;
    try std.testing.expectEqual(@as(?usize, null), try engine.simplex(3, false));
    try std.testing.expectEqual(@as(usize, 0), w.basis[1]);
    try std.testing.expectEqual(@as(usize, 1), engine.iterations);
    try std.testing.expectApproxEqAbs(@as(f64, 2e-15), engine.cell(1, 5).*, 1e-20);
    try std.testing.expect(@abs(engine.cell(0, 5).*) < 1e-20);
}
