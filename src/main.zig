const std = @import("std");
const RootFinding = @import("rootfinding.zig").rootfinding;
const NumericalInteg = @import("numericalinteg.zig").numericalinteg;

// Example functions for root finding
fn f1(x: f32) f32 {
    return x * x - 4.0; // Root at x = 2.0
}

fn df1(x: f32) f32 {
    return 2.0 * x;
}

fn g1(x: f32) f32 {
    return x - (x * x - 4.0) / 10.0; // Fixed point form of f1
}

fn f2(x: f32) f32 {
    return x * x * x - x - 1.0; // Root near x = 1.32
}

fn df2(x: f32) f32 {
    return 3.0 * x * x - 1.0;
}

fn g2(x: f32) f32 {
    return std.math.pow(f32, x + 1.0, 1.0 / 3.0); // Fixed point form of f2
}

// Example function for integration: f(x) = x^2
fn integ(x: f32) f32 {
    return x * x;
}

pub fn main() !void {
    const tolerance: f32 = 0.0001;
    const max_iter: usize = 100;

    // Root finding examples
    const root_bisection = try RootFinding.bisection(f1, 0.0, 3.0, tolerance, max_iter);
    std.debug.print("Bisection root: {d:.6}\n", .{root_bisection});

    const root_falsepos = try RootFinding.falsePosition(f1, 0.0, 3.0, tolerance, max_iter);
    std.debug.print("False Position root: {d:.6}\n", .{root_falsepos});

    const root_fixedpoint = try RootFinding.fixedPoint(g1, 1.0, tolerance, max_iter);
    std.debug.print("Fixed Point root: {d:.6}\n", .{root_fixedpoint});

    const root_newton = try RootFinding.newtonRaphson(f2, df2, 1.0, tolerance, max_iter);
    std.debug.print("Newton-Raphson root: {d:.6}\n", .{root_newton});

    const root_brent = try RootFinding.brent(f2, 1.0, 2.0, tolerance, max_iter);
    std.debug.print("Brent's method root: {d:.6}\n", .{root_brent});

    // Numerical integration examples
    const a: f32 = 0.0;
    const b: f32 = 1.0;
    const n_even: usize = 10; // For simpson13 (must be even)
    const n_mult3: usize = 12; // For simpson38 (must be multiple of 3)
    const n: usize = 10; // For trapezoidal

    const trap = try NumericalInteg.trapezoidal(integ, a, b, n);
    std.debug.print("Trapezoidal Rule: {d:.6}\n", .{trap});

    const simp13 = try NumericalInteg.simpson13(integ, a, b, n_even);
    std.debug.print("Simpson's 1/3 Rule: {d:.6}\n", .{simp13});

    const simp38 = try NumericalInteg.simpson38(integ, a, b, n_mult3);
    std.debug.print("Simpson's 3/8 Rule: {d:.6}\n", .{simp38});

    const gauss = try NumericalInteg.gaussianQuadrature(integ, a, b);
    std.debug.print("Gaussian Quadrature (2-point): {d:.6}\n", .{gauss});
}
