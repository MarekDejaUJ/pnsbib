const std = @import("std");

pub const Bound = struct {
    lower: f64,
    upper: f64,
};

/// Closed-form Li-Pearl bounds for Theorems 4-7.
pub fn single(
    kind: u8,
    aji: f64,
    oji: f64,
    vi: f64,
    vk: f64,
    uj: f64,
    up: f64,
    opk: f64,
    ojk: f64,
    other_ok: []const f64,
) Bound {
    switch (kind) {
        4 => return .{
            .lower = @max(oji, aji + vi - 1),
            .upper = @min(aji, vi),
        },
        5 => {
            var partition: f64 = 0;
            for (other_ok) |value| {
                partition += @max(0, aji + value - 1 + uj - oji);
            }
            return .{
                .lower = @max(0, @max(aji + vk - 1, partition)),
                .upper = @min(aji - oji, vk - ojk),
            };
        },
        6 => return .{
            .lower = @max(0, aji - oji - 1 + uj + up),
            .upper = @min(aji - oji, up),
        },
        7 => return .{
            .lower = @max(0, aji + opk - 1 + uj - oji),
            .upper = @min(aji - oji, opk),
        },
        else => unreachable,
    }
}

test "preservation retains the factual lower bound" {
    const result = single(4, 0.4, 0.2, 0.5, 0, 0.6, 0, 0, 0, &.{});
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), result.lower, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), result.upper, 1e-12);
}

test "substitute respects available mass in the factual arm" {
    const result = single(6, 0.7, 0.2, 0, 0, 0.4, 0.3, 0, 0, &.{});
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), result.lower, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), result.upper, 1e-12);
}
