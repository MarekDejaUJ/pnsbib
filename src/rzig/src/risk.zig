const std = @import("std");

pub const Error = error{ InvalidDimensions, InvalidInput, ResourceLimit, Singular, NoConvergence, NumericalFailure };
pub const Result = struct { iterations: usize, objective: f64, gradient: f64 };

pub fn workspaceSize(n: usize, p: usize, m: usize) Error!usize {
    if (n == 0 or p == 0 or m == 0) return error.InvalidDimensions;
    if (p > 256 or n > 8_000_000 / p or m > 8_000_000 / p) return error.ResourceLimit;
    return p * p + 5 * p;
}

fn sigmoid(z: f64) f64 {
    if (z >= 0) return 1 / (1 + @exp(-z));
    const e = @exp(z);
    return e / (1 + e);
}

fn dot(n: usize, p: usize, x: []const f64, i: usize, beta: []const f64) f64 {
    var z: f64 = 0;
    for (0..p) |j| z += x[i + n * j] * beta[j];
    return z;
}

fn objective(n: usize, p: usize, x: []const f64, y: []const f64, w: []const f64, penalty: []const f64, lambda: f64, beta: []const f64) f64 {
    var value: f64 = 0;
    for (0..n) |i| {
        const z = dot(n, p, x, i, beta);
        value += w[i] * (@max(z, 0) - y[i] * z + std.math.log1p(@exp(-@abs(z))));
    }
    for (0..p) |j| value += 0.5 * lambda * penalty[j] * beta[j] * beta[j];
    return value;
}

/// Strictly convex penalized binomial fit. All matrices are column-major.
/// Caller owns scratch; no allocation, mutable global state or R fallback.
pub fn fit(n: usize, p: usize, m: usize, x: []const f64, y: []const f64, w: []const f64, xp: []const f64, penalty: []const f64, lambda: f64, max_iterations: usize, tolerance: f64, beta: []f64, pred: []f64, scratch: []f64) Error!Result {
    const size = try workspaceSize(n, p, m);
    if (x.len != n * p or y.len != n or w.len != n or xp.len != m * p or penalty.len != p or beta.len != p or pred.len != m or scratch.len < size) return error.InvalidDimensions;
    if (!std.math.isFinite(lambda) or lambda <= 0 or !std.math.isFinite(tolerance) or tolerance <= 0 or tolerance > 1e-3 or max_iterations == 0) return error.InvalidInput;
    for (x) |v| if (!std.math.isFinite(v)) {
        return error.InvalidInput;
    };
    for (xp) |v| if (!std.math.isFinite(v)) {
        return error.InvalidInput;
    };
    var sumw: f64 = 0;
    for (0..n) |i| {
        if (!std.math.isFinite(y[i]) or y[i] < 0 or y[i] > 1 or !std.math.isFinite(w[i]) or w[i] < 0) return error.InvalidInput;
        sumw += w[i];
    }
    if (!std.math.isFinite(sumw) or sumw <= 0) return error.InvalidInput;
    for (penalty) |v| if (!std.math.isFinite(v) or v <= 0) {
        return error.InvalidInput;
    };
    const h = scratch[0 .. p * p];
    const grad = scratch[p * p .. p * p + p];
    const step = scratch[p * p + p .. p * p + 2 * p];
    const candidate = scratch[p * p + 2 * p .. p * p + 3 * p];
    @memset(beta, 0);
    var value = objective(n, p, x, y, w, penalty, lambda, beta);
    for (0..max_iterations) |iter| {
        @memset(h, 0);
        for (0..p) |j| {
            grad[j] = lambda * penalty[j] * beta[j];
            h[j * p + j] = lambda * penalty[j];
        }
        for (0..n) |i| {
            const q = sigmoid(dot(n, p, x, i, beta));
            for (0..p) |j| {
                grad[j] += w[i] * (q - y[i]) * x[i + n * j];
                for (0..j + 1) |k| h[j * p + k] += w[i] * q * (1 - q) * x[i + n * j] * x[i + n * k];
            }
        }
        var norm: f64 = 0;
        for (grad) |v| norm = @max(norm, @abs(v));
        norm /= @max(1, sumw);
        if (!std.math.isFinite(norm) or !std.math.isFinite(value)) return error.NumericalFailure;
        if (norm <= tolerance) {
            for (0..m) |i| pred[i] = sigmoid(dot(m, p, xp, i, beta));
            return .{ .iterations = iter, .objective = value, .gradient = norm };
        }
        // Cholesky factor of the penalized Hessian, then two triangular solves.
        for (0..p) |j| {
            for (0..j + 1) |k| {
                var v = h[j * p + k];
                for (0..k) |l| v -= h[j * p + l] * h[k * p + l];
                if (j == k) {
                    if (!std.math.isFinite(v) or v <= 0) return error.Singular;
                    h[j * p + k] = @sqrt(v);
                } else h[j * p + k] = v / h[k * p + k];
            }
            var v = grad[j];
            for (0..j) |k| v -= h[j * p + k] * step[k];
            step[j] = v / h[j * p + j];
        }
        for (0..p) |back| {
            const j = p - 1 - back;
            var v = step[j];
            for (j + 1..p) |k| v -= h[k * p + j] * step[k];
            step[j] = v / h[j * p + j];
        }
        var descent: f64 = 0;
        for (0..p) |j| descent += grad[j] * step[j];
        var accepted = false;
        var alpha: f64 = 1;
        for (0..50) |_| {
            for (0..p) |j| candidate[j] = beta[j] - alpha * step[j];
            const next = objective(n, p, x, y, w, penalty, lambda, candidate);
            if (std.math.isFinite(next) and next <= value - 1e-4 * alpha * descent + 1e-12 * @max(1, value)) {
                @memcpy(beta, candidate);
                value = next;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) return error.NumericalFailure;
    }
    return error.NoConvergence;
}

pub fn standardize(n: usize, t: usize, pred: []const f64, weights: []const f64, cumulative: bool, out: []f64) Error!void {
    if (n == 0 or t == 0 or n > 8_000_000 or t > 8_000_000 / n) return error.InvalidDimensions;
    if (pred.len != n * t or weights.len != n or out.len != t) return error.InvalidDimensions;
    var sum: f64 = 0;
    for (weights) |v| {
        if (!std.math.isFinite(v) or v < 0) return error.InvalidInput;
        sum += v;
    }
    if (!std.math.isFinite(sum) or sum <= 0) return error.InvalidInput;
    for (pred) |v| if (!std.math.isFinite(v) or v < 0 or v > 1) {
        return error.InvalidInput;
    };
    @memset(out, 0);
    for (0..n) |i| {
        var survival: f64 = 1;
        for (0..t) |h| {
            const q = pred[i + n * h];
            survival *= 1 - q;
            out[h] += weights[i] / sum * (if (cumulative) 1 - survival else q);
        }
    }
}

test "balanced intercept has zero coefficient and half risk" {
    var beta: [1]f64 = undefined;
    var pred: [2]f64 = undefined;
    var scratch: [6]f64 = undefined;
    const r = try fit(2, 1, 2, &.{ 1, 1 }, &.{ 0, 1 }, &.{ 1, 1 }, &.{ 1, 1 }, &.{0.1}, 1, 100, 1e-10, &beta, &pred, &scratch);
    try std.testing.expectApproxEqAbs(@as(f64, 0), beta[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), pred[0], 1e-12);
    try std.testing.expect(r.gradient <= 1e-10);
}

test "all successes remain finite with positive penalty" {
    var beta: [1]f64 = undefined;
    var pred: [1]f64 = undefined;
    var scratch: [6]f64 = undefined;
    _ = try fit(2, 1, 1, &.{ 1, 1 }, &.{ 1, 1 }, &.{ 1, 1 }, &.{1}, &.{0.1}, 1, 100, 1e-10, &beta, &pred, &scratch);
    try std.testing.expect(pred[0] > 0.8 and pred[0] < 1);
    try std.testing.expectApproxEqAbs(2 * (pred[0] - 1) + 0.1 * beta[0], 0, 1e-8);
}

test "invalid outcomes fail explicitly" {
    var beta: [1]f64 = undefined;
    var pred: [1]f64 = undefined;
    var scratch: [6]f64 = undefined;
    try std.testing.expectError(error.InvalidInput, fit(1, 1, 1, &.{1}, &.{2}, &.{1}, &.{1}, &.{1}, 1, 100, 1e-10, &beta, &pred, &scratch));
}

test "cumulative standardization uses survival products" {
    var out: [2]f64 = undefined;
    try standardize(2, 2, &.{ 0.2, 0.4, 0.5, 0.2 }, &.{ 1, 3 }, true, &out);
    try std.testing.expectApproxEqAbs(@as(f64, 0.35), out[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.54), out[1], 1e-12);
}
