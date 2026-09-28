const std = @import("std");

pub const Error = error{ InvalidDimensions, InvalidInput, InvalidQuery, VariableLimit, WorkspaceLimit, Interrupted };
pub const Poll = *const fn () error{Interrupted}!void;
pub const max_cells: usize = 16_777_216;
pub const Shape = struct { variables: usize, rows: usize, types: usize, cells: usize };

fn mul(a: usize, b: usize) Error!usize {
    return std.math.mul(usize, a, b) catch error.WorkspaceLimit;
}

fn add(a: usize, b: usize) Error!usize {
    return std.math.add(usize, a, b) catch error.WorkspaceLimit;
}

pub fn shape(m: usize, n: usize, max_variables: usize) Error!Shape {
    if (m < 2 or n < 2 or max_variables == 0) return error.InvalidDimensions;
    if (m > std.math.maxInt(i32) or n > std.math.maxInt(i32)) return error.WorkspaceLimit;
    var types: usize = 1;
    for (0..m) |_| {
        if (types > max_variables / n) return error.VariableLimit;
        types *= n;
    }
    if (types > max_variables / m) return error.VariableLimit;
    const variables = try mul(m, types);
    const rows = try add(1, try mul(2, try mul(m, n)));
    const cells = try add(try add(3, rows), try mul(variables, try add(try add(m, 2), rows)));
    if (cells > max_cells) return error.WorkspaceLimit;
    return .{ .variables = variables, .rows = rows, .types = types, .cells = cells };
}

/// Output: variables, rows, types, response (column-major), factual X, factual
/// Y, A (column-major), b. Margins are copied unchanged; R validates their
/// probability/causal interpretation. All caller-owned memory is setup-only.
pub fn build(m: usize, n: usize, observed: []const f64, intervention: []const f64, max_variables: usize, out: []f64, poll: ?Poll) Error!void {
    const s = try shape(m, n, max_variables);
    if (observed.len != m * n or intervention.len != m * n or out.len != s.cells) return error.InvalidDimensions;
    for (observed) |v| if (!std.math.isFinite(v)) return error.InvalidInput;
    for (intervention) |v| if (!std.math.isFinite(v)) return error.InvalidInput;
    if (poll) |callback| try callback();
    out[0] = @floatFromInt(s.variables);
    out[1] = @floatFromInt(s.rows);
    out[2] = @floatFromInt(s.types);
    const response = out[3 .. 3 + s.variables * m];
    const fx = out[3 + s.variables * m .. 3 + s.variables * (m + 1)];
    const fy = out[3 + s.variables * (m + 1) .. 3 + s.variables * (m + 2)];
    const a = out[3 + s.variables * (m + 2) .. s.cells - s.rows];
    const b = out[s.cells - s.rows ..];
    @memset(a, 0);
    b[0] = 1;
    for (0..m) |x| for (0..n) |y| {
        b[1 + x * n + y] = observed[x + y * m];
        b[1 + m * n + x * n + y] = intervention[x + y * m];
    };
    for (0..s.variables) |v| {
        if (v % 256 == 0) if (poll) |callback| try callback();
        const factual = v / s.types;
        fx[v] = @floatFromInt(factual + 1);
        var code = v % s.types;
        var factual_y: usize = 0;
        a[v * s.rows] = 1;
        for (0..m) |x| {
            const y = code % n;
            code /= n;
            response[v + x * s.variables] = @floatFromInt(y + 1);
            if (x == factual) factual_y = y;
            a[1 + m * n + x * n + y + v * s.rows] = 1;
        }
        fy[v] = @floatFromInt(factual_y + 1);
        a[1 + factual * n + factual_y + v * s.rows] = 1;
    }
    if (poll) |callback| try callback();
}

pub fn eventSize(variables: usize, m: usize) Error!usize {
    if (variables == 0 or m < 2) return error.InvalidDimensions;
    if (m > std.math.maxInt(i32) or try mul(variables, m) > max_cells) return error.WorkspaceLimit;
    return m;
}

fn level(value: i32, maximum: usize) bool {
    return value >= 1 and @as(usize, @intCast(value)) <= maximum;
}

/// Positive integer category codes need not be contiguous. Factual Y must
/// equal the potential response selected by factual X. Zero denotes an absent
/// factual condition; the empty counterfactual conjunction is supported.
pub fn eventMask(variables: usize, m: usize, response: []const i32, fx: []const i32, fy: []const i32, x: []const i32, y: []const i32, ox: i32, oy: i32, out: []f64, selected: []bool, poll: ?Poll) Error!void {
    _ = try eventSize(variables, m);
    if (response.len != variables * m or fx.len != variables or fy.len != variables or out.len != variables or selected.len != m) return error.InvalidDimensions;
    if (x.len != y.len or x.len > m or ox < 0 or (ox != 0 and !level(ox, m)) or oy < 0) return error.InvalidQuery;
    if (poll) |callback| try callback();
    @memset(selected, false);
    for (x, y, 0..) |xx, yy, i| {
        if (i % 1024 == 0) if (poll) |callback| try callback();
        if (!level(xx, m) or yy < 1) return error.InvalidQuery;
        const index: usize = @intCast(xx - 1);
        if (selected[index]) return error.InvalidQuery;
        selected[index] = true;
    }
    for (response, 0..) |value, i| {
        if (i % 4096 == 0) if (poll) |callback| try callback();
        if (value < 1) return error.InvalidInput;
    }
    for (0..variables) |v| {
        if (v % 256 == 0) if (poll) |callback| try callback();
        if (!level(fx[v], m) or fy[v] < 1) return error.InvalidInput;
        if (fy[v] != response[v + @as(usize, @intCast(fx[v] - 1)) * variables]) return error.InvalidInput;
        var keep = (ox == 0 or fx[v] == ox) and (oy == 0 or fy[v] == oy);
        for (x, y) |xx, yy| keep = keep and response[v + @as(usize, @intCast(xx - 1)) * variables] == yy;
        out[v] = if (keep) 1 else 0;
    }
    if (poll) |callback| try callback();
}

test "binary static cells preserve factual-block and potential-digit order" {
    const s = try shape(2, 2, 10000);
    try std.testing.expectEqual(@as(usize, 8), s.variables);
    try std.testing.expectEqual(@as(usize, 9), s.rows);
    var out: [116]f64 = undefined;
    try build(2, 2, &.{ 0.3, 0.1, 0.2, 0.4 }, &.{ 0.6, 0.2, 0.4, 0.8 }, 10000, &out, null);
    try std.testing.expectEqualSlices(f64, &.{ 8, 9, 4 }, out[0..3]);
    try std.testing.expectEqualSlices(f64, &.{ 1, 2, 1, 2, 1, 2, 1, 2 }, out[3..11]);
    try std.testing.expectEqualSlices(f64, &.{ 1, 1, 2, 2, 1, 1, 2, 2 }, out[11..19]);
    try std.testing.expectEqualSlices(f64, &.{ 1, 1, 1, 1, 2, 2, 2, 2 }, out[19..27]);
    try std.testing.expectEqualSlices(f64, &.{ 1, 2, 1, 2, 1, 1, 2, 2 }, out[27..35]);
    try std.testing.expectEqualSlices(f64, &.{ 1, 1, 0, 0, 0, 1, 0, 1, 0 }, out[35..44]);
    try std.testing.expectEqualSlices(f64, &.{ 1, 0.3, 0.2, 0.1, 0.4, 0.6, 0.4, 0.2, 0.8 }, out[107..116]);
}

test "static shape refuses oversized work before multiplication or allocation" {
    try std.testing.expectError(error.InvalidDimensions, shape(1, 2, 10));
    try std.testing.expectError(error.VariableLimit, shape(2, 2, 7));
    try std.testing.expectError(error.VariableLimit, shape(1000, 2, 10000));
    try std.testing.expectError(error.WorkspaceLimit, shape(2, 1000, 10000000));
    try std.testing.expectError(error.WorkspaceLimit, shape(std.math.maxInt(usize), 2, std.math.maxInt(usize)));
}

test "static construction rejects nonfinite margins and preserves borrowed values" {
    var out: [116]f64 = undefined;
    const o = [_]f64{ 0.3, 0.1, 0.2, 0.4 };
    try std.testing.expectError(error.InvalidDimensions, build(2, 2, &o, &.{0.5}, 10000, &out, null));
    try std.testing.expectError(error.InvalidInput, build(2, 2, &o, &.{ 0, 1, std.math.nan(f64), 1 }, 10000, &out, null));
    try std.testing.expectEqualSlices(f64, &.{ 0.3, 0.1, 0.2, 0.4 }, &o);
}

test "categorical event masks include factual conditions and empty queries" {
    var out: [2]f64 = undefined;
    var scratch: [2]bool = undefined;
    try eventMask(2, 2, &.{ 1, 2, 2, 1 }, &.{ 1, 2 }, &.{ 1, 1 }, &.{1}, &.{2}, 0, 0, &out, &scratch, null);
    try std.testing.expectEqualSlices(f64, &.{ 0, 1 }, &out);
    try eventMask(2, 2, &.{ 1, 2, 2, 1 }, &.{ 1, 2 }, &.{ 1, 1 }, &.{}, &.{}, 0, 1, &out, &scratch, null);
    try std.testing.expectEqualSlices(f64, &.{ 1, 1 }, &out);
    try eventMask(2, 2, &.{ 1, 2, 2, 1 }, &.{ 1, 2 }, &.{ 1, 1 }, &.{1}, &.{1}, 2, 2, &out, &scratch, null);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0 }, &out);
    try std.testing.expectError(error.InvalidQuery, eventMask(2, 2, &.{ 1, 2, 2, 1 }, &.{ 1, 2 }, &.{ 1, 1 }, &.{ 1, 1 }, &.{ 1, 1 }, 0, 0, &out, &scratch, null));
    try std.testing.expectError(error.InvalidInput, eventMask(2, 2, &.{ 1, 2, 2, 1 }, &.{ 1, 3 }, &.{ 1, 1 }, &.{}, &.{}, 0, 0, &out, &scratch, null));
}

test "static numerical kernels propagate interrupts" {
    const Interrupt = struct {
        fn poll() error{Interrupted}!void {
            return error.Interrupted;
        }
    };
    var out: [116]f64 = undefined;
    try std.testing.expectError(error.Interrupted, build(2, 2, &.{ 0.3, 0.1, 0.2, 0.4 }, &.{ 0.6, 0.2, 0.4, 0.8 }, 10000, &out, Interrupt.poll));
    var mask: [2]f64 = undefined;
    var scratch: [2]bool = undefined;
    try std.testing.expectError(error.Interrupted, eventMask(2, 2, &.{ 1, 2, 2, 1 }, &.{ 1, 2 }, &.{ 1, 1 }, &.{}, &.{}, 0, 0, &mask, &scratch, Interrupt.poll));
}
