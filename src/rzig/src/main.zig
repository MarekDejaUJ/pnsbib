const std = @import("std");
const builtin = @import("builtin");
const rzig = @import("rzig");
const bounds = @import("bounds.zig");
const paper = @import("paper.zig");
const shu = @import("shu.zig");
const lp = @import("lp.zig");
const static = @import("static.zig");
const risk = @import("risk.zig");

/// Internal penalized binary-risk fit; column-major numeric inputs.
pub fn poc_risk_fit_zig(ctx: *rzig.Ctx, n: usize, p: usize, m: usize, x: []const f64, y: []const f64, weights: []const f64, target: []const f64, penalty: []const f64, lambda: f64, max_iterations: usize, tolerance: f64) rzig.Error![]f64 {
    try rzig.checkInterrupt();
    const size = risk.workspaceSize(n, p, m) catch |err| return rzig.raise("Risk dimensions: {s}", .{@errorName(err)});
    const output = try ctx.alloc(f64, 3 + p + m);
    const scratch = try ctx.alloc(f64, size);
    const result = risk.fit(n, p, m, x, y, weights, target, penalty, lambda, max_iterations, tolerance, output[3 .. 3 + p], output[3 + p ..], scratch) catch |err| return rzig.raise("Risk fit: {s}", .{@errorName(err)});
    output[0] = result.objective;
    output[1] = @floatFromInt(result.iterations);
    output[2] = result.gradient;
    try rzig.checkInterrupt();
    return output;
}

/// Internal target standardization, returning risks and individual paths.
pub fn poc_risk_standardize_zig(ctx: *rzig.Ctx, n: usize, t: usize, pred: []const f64, weights: []const f64, cumulative: bool) rzig.Error![]f64 {
    if (n == 0 or t == 0 or n > 8_000_000 or t > 8_000_000 / n) return rzig.raise("Risk dimensions: InvalidDimensions", .{});
    const output = try ctx.alloc(f64, t + n * t);
    risk.standardize(n, t, pred, weights, cumulative, output[0..t]) catch |err| return rzig.raise("Risk standardization: {s}", .{@errorName(err)});
    for (0..n) |i| {
        var survival: f64 = 1;
        for (0..t) |h| {
            const q = pred[i + n * h];
            survival *= 1 - q;
            output[t + i + n * h] = if (cumulative) 1 - survival else q;
        }
    }
    return output;
}

/// Route ReleaseSafe failures through R while retaining Zig test behavior.
pub const panic = if (builtin.is_test)
    std.debug.FullPanic(std.debug.defaultPanic)
else
    rzig.Panic;

/// Compute one single-term Li-Pearl probability-of-causation interval.
/// @param kind Theorem number, 4 through 7.
/// @param aji Interventional probability for the target treatment and outcome.
/// @param oji Observed joint probability for the target treatment and outcome.
/// @param vi Observed marginal probability of the target outcome.
/// @param vk Observed marginal probability of the factual outcome.
/// @param uj Observed probability of the target treatment.
/// @param up Observed probability of the factual treatment.
/// @param opk Observed joint probability of factual treatment and outcome.
/// @param ojk Observed joint probability of target treatment and factual outcome.
/// @param other_ok Observed factual-outcome probabilities for other treatments.
/// @return Lower and upper bounds.
/// @export
pub fn poc_single_zig(
    ctx: *rzig.Ctx,
    kind: i32,
    aji: f64,
    oji: f64,
    vi: f64,
    vk: f64,
    uj: f64,
    up: f64,
    opk: f64,
    ojk: f64,
    other_ok: []const f64,
) rzig.Error![]f64 {
    if (kind < 4 or kind > 7) return rzig.raise("kind must be 4, 5, 6 or 7", .{});
    const inputs = [_]f64{ aji, oji, vi, vk, uj, up, opk, ojk };
    for (inputs) |value| {
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return rzig.raise("probabilities must be finite and in [0, 1]", .{});
    }
    for (other_ok) |value| {
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return rzig.raise("probabilities must be finite and in [0, 1]", .{});
    }
    const result = bounds.single(@intCast(kind), aji, oji, vi, vk, uj, up, opk, ojk, other_ok);
    const output = try ctx.alloc(f64, 2);
    output[0] = result.lower;
    output[1] = result.upper;
    return output;
}

fn pollPaper() error{Interrupted}!void {
    rzig.checkInterrupt() catch return error.Interrupted;
}

/// Evaluate the complete printed Li-Pearl analytical recursion in Zig.
/// @param m Number of treatment levels, at least two.
/// @param n Number of outcome levels, at least two.
/// @param observed Joint probabilities in column-major matrix order.
/// @param intervention Intervention probabilities in column-major order.
/// @param x Distinct one-based counterfactual treatment indices.
/// @param y Corresponding one-based outcome indices.
/// @param observed_x One-based factual treatment index, or zero for none.
/// @param observed_y One-based factual outcome index, or zero for none.
/// @param tolerance Positive probability-consistency tolerance.
/// @return Unconditional lower and upper bounds. The memo is limited to 1048576 states.
/// @export
pub fn poc_paper_zig(
    ctx: *rzig.Ctx,
    m: usize,
    n: usize,
    observed: []const f64,
    intervention: []const f64,
    x: []const i32,
    y: []const i32,
    observed_x: i32,
    observed_y: i32,
    tolerance: f64,
) rzig.Error![]f64 {
    const size = paper.workspaceSize(m, n, x.len) catch |err|
        return rzig.raise("Li-Pearl native input: {s}", .{@errorName(err)});
    const px = try ctx.alloc(f64, m);
    const py = try ctx.alloc(f64, n);
    const memo = try ctx.alloc(paper.Bound, size);
    const answer = paper.evaluateInterruptible(m, n, observed, intervention, x, y, observed_x, observed_y, tolerance, px, py, memo, pollPaper) catch |err| {
        if (err == error.Interrupted) return error.RZigError;
        return rzig.raise("Li-Pearl native input: {s}", .{@errorName(err)});
    };
    const output = try ctx.alloc(f64, 2);
    output[0] = answer.lower;
    output[1] = answer.upper;
    return output;
}

/// Evaluate the Shu-Wang-Li closed form and the labelled repeated-outcome extension.
/// @param m Number of treatment levels, at least two.
/// @param n Number of outcome levels, at least two.
/// @param observed Joint probabilities in column-major matrix order.
/// @param intervention Intervention probabilities in column-major order.
/// @param x Distinct one-based counterfactual treatment indices.
/// @param y Corresponding one-based outcome indices.
/// @param observed_x One-based factual treatment index, or zero for none.
/// @param observed_y One-based factual outcome index, or zero for none.
/// @param tolerance Positive probability-consistency tolerance.
/// @return Lower and upper bounds, family code and consistency-reduction flag.
/// Family codes 0-6 denote factual margin, contradiction, theorems 2-5, and the
/// repeated-outcome extension. The complete input table has at most 1048576 cells.
/// @export
pub fn poc_shu_zig(
    ctx: *rzig.Ctx,
    m: usize,
    n: usize,
    observed: []const f64,
    intervention: []const f64,
    x: []const i32,
    y: []const i32,
    observed_x: i32,
    observed_y: i32,
    tolerance: f64,
) rzig.Error![]f64 {
    try rzig.checkInterrupt();
    shu.checkSize(m, n) catch |err|
        return rzig.raise("Shu native input: {s}", .{@errorName(err)});
    if (x.len > m or y.len != x.len) return rzig.raise("Shu native input: InvalidQuery", .{});
    const px = try ctx.alloc(f64, m);
    const py = try ctx.alloc(f64, n);
    const selected = try ctx.alloc(bool, m);
    const differences = try ctx.alloc(f64, x.len);
    const result = shu.evaluate(m, n, observed, intervention, x, y, observed_x, observed_y, tolerance, px, py, selected, differences) catch |err|
        return rzig.raise("Shu native input: {s}", .{@errorName(err)});
    try rzig.checkInterrupt();
    const output = try ctx.alloc(f64, 4);
    output[0] = result.bounds.lower;
    output[1] = result.bounds.upper;
    output[2] = @floatFromInt(@intFromEnum(result.family));
    output[3] = @floatFromInt(@intFromBool(result.after_consistency));
    return output;
}

/// Solve a continuous nonnegative-variable LP using the native two-phase solver.
/// @param m Number of constraint rows, possibly zero.
/// @param n Number of variables, at least one.
/// @param a Constraint matrix in column-major order.
/// @param b Constraint right-hand sides.
/// @param directions Integer directions: -1 for less-than, 0 equality, 1 greater-than.
/// @param objective Objective coefficients.
/// @param maximize Logical maximization flag.
/// @param max_iterations Positive pivot limit, normally 200000.
/// @return Numeric vector: status (0 optimal, 2 infeasible, 3 unbounded), objective,
/// iterations, primal violation, dual violation, duality gap, then n primal values,
/// m dual values and n ray values. Infeasible dual values are a Farkas certificate;
/// an unbounded result contains a feasible point and improving ray. Nonoptimal
/// objective values are NaN. Numerical failures and resource limits raise errors.
/// Pivot and certificate tolerances are 1e-10 and 1e-8; tableau plus lexicographic
/// scratch is capped at 16777216 cells. Pivoting uses an improving coefficient
/// and symbolic lexicographic row ordering without changing constraint values.
/// @export
pub fn poc_lp_zig(
    ctx: *rzig.Ctx,
    m: usize,
    n: usize,
    a: []const f64,
    b: []const f64,
    directions: []const i32,
    objective: []const f64,
    maximize: bool,
    max_iterations: usize,
) rzig.Error![]f64 {
    try rzig.checkInterrupt();
    _ = lp.workspaceSize(m, n) catch |err|
        return rzig.raise("LP native input: {s}", .{@errorName(err)});
    if (a.len != m * n or b.len != m or directions.len != m or objective.len != n)
        return rzig.raise("LP native input: InvalidDimensions", .{});
    var workspace = lp.Workspace.initAlloc(ctx.allocator(), m, n) catch |err|
        return rzig.raise("LP native setup: {s}", .{@errorName(err)});
    const output = try ctx.alloc(f64, 6 + 2 * n + m);
    const primal = output[6 .. 6 + n];
    const dual = output[6 + n .. 6 + n + m];
    const ray = output[6 + n + m ..];
    const result = lp.solve(m, n, a, b, directions, objective, maximize, .{ .max_iterations = max_iterations }, primal, dual, ray, &workspace, pollPaper) catch |err| {
        if (err == error.Interrupted) return error.RZigError;
        return rzig.raise("LP native solve: {s}", .{@errorName(err)});
    };
    output[0] = switch (result.status) {
        .optimal => 0,
        .infeasible => 2,
        .unbounded => 3,
    };
    output[1] = result.objective;
    output[2] = @floatFromInt(result.iterations);
    output[3] = result.max_primal_violation;
    output[4] = result.max_dual_violation;
    output[5] = result.duality_gap;
    try rzig.checkInterrupt();
    return output;
}

/// Construct the complete static response-type system in native code.
/// @param m Number of treatment levels, at least two.
/// @param n Number of outcome levels, at least two.
/// @param observed Observational margins in column-major order, copied unchanged.
/// @param intervention Interventional margins in column-major order, copied unchanged.
/// @param max_variables Positive response-variable cap.
/// @return Numeric vector: variable count V, row count K, response-type count,
/// then V*m response codes, V factual-X codes, V factual-Y codes, K*V constraint
/// coefficients and K right-hand sides. Matrices are column-major. Output is
/// capped at 16777216 doubles. R retains probability/label validation.
/// @export
pub fn poc_static_system_zig(ctx: *rzig.Ctx, m: usize, n: usize, observed: []const f64, intervention: []const f64, max_variables: usize) rzig.Error![]f64 {
    try rzig.checkInterrupt();
    const s = static.shape(m, n, max_variables) catch |err|
        return rzig.raise("Static native shape: {s}", .{@errorName(err)});
    if (observed.len != m * n or intervention.len != m * n)
        return rzig.raise("Static native input: InvalidDimensions", .{});
    const output = try ctx.alloc(f64, s.cells);
    static.build(m, n, observed, intervention, max_variables, output, pollPaper) catch |err| {
        if (err == error.Interrupted) return error.RZigError;
        return rzig.raise("Static native construction: {s}", .{@errorName(err)});
    };
    return output;
}

/// Evaluate a categorical response-type conjunction and factual condition.
/// @param variables Number of response-distribution cells.
/// @param m Number of treatment columns, at least two.
/// @param response Positive integer potential-outcome codes, column-major.
/// @param factual_x One-based factual treatment indices.
/// @param factual_y Factual outcome codes consistent with the response matrix.
/// @param x Distinct one-based counterfactual treatment indices; may be empty.
/// @param y Corresponding positive outcome codes.
/// @param observed_x One-based factual treatment or zero for no condition.
/// @param observed_y Positive factual outcome code or zero for no condition.
/// @return Numeric zero/one event mask. At most 16777216 response cells.
/// @export
pub fn poc_static_event_zig(ctx: *rzig.Ctx, variables: usize, m: usize, response: []const i32, factual_x: []const i32, factual_y: []const i32, x: []const i32, y: []const i32, observed_x: i32, observed_y: i32) rzig.Error![]f64 {
    try rzig.checkInterrupt();
    const size = static.eventSize(variables, m) catch |err|
        return rzig.raise("Event native shape: {s}", .{@errorName(err)});
    if (response.len != variables * m or factual_x.len != variables or factual_y.len != variables)
        return rzig.raise("Event native input: InvalidDimensions", .{});
    const selected = try ctx.alloc(bool, size);
    const output = try ctx.alloc(f64, variables);
    static.eventMask(variables, m, response, factual_x, factual_y, x, y, observed_x, observed_y, output, selected, pollPaper) catch |err| {
        if (err == error.Interrupted) return error.RZigError;
        return rzig.raise("Event native input: {s}", .{@errorName(err)});
    };
    return output;
}

comptime {
    rzig.registerModule(@This());
}
