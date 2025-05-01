const std = @import("std");
const GaussSeidel = @import("gaussseidel.zig").gaussseidel;
const RootFinding = @import("rootfinding.zig").rootfinding;
const NumericalInteg = @import("numericalinteg.zig").numericalinteg;
const LPP = @import("lpp.zig").lpp;

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

// Example functions for integration
fn integ1(x: f32) f32 {
    return x * x; // ∫x² dx = x³/3
}

fn integ2(x: f32) f32 {
    return std.math.sin(x); // ∫sin(x) dx = -cos(x)
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Parse command line arguments
    var args = try std.process.argsWithAllocator(allocator);
    defer args.deinit();
    _ = args.skip(); // Skip program name

    var show_iterations = false;
    var method: []const u8 = "gaussseidel";

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--show-iterations")) {
            show_iterations = true;
        } else if (std.mem.eql(u8, arg, "--method")) {
            if (args.next()) |m| {
                method = m;
            } else {
                std.debug.print("Error: --method requires a value\n", .{});
                return;
            }
        } else if (std.mem.eql(u8, arg, "--help")) {
            std.debug.print("Usage: main [options]\n", .{});
            std.debug.print("Options:\n", .{});
            std.debug.print("  --show-iterations  Show iteration table\n", .{});
            std.debug.print("  --method METHOD    Choose method:\n", .{});
            std.debug.print("                     - gaussseidel\n", .{});
            std.debug.print("                     - bisection, falsepos, fixedpoint, newton, brent\n", .{});
            std.debug.print("                     - trapezoidal, simpson, romberg, gaussian\n", .{});
            std.debug.print("                     - lpp\n", .{});
            std.debug.print("  --help            Show this help message\n", .{});
            return;
        }
    }

    if (std.mem.eql(u8, method, "gaussseidel")) {
        // Run Gauss-Seidel method
        const n: usize = 5;
        var matrix_data = try allocator.alloc([]f32, n);
        errdefer {
            for (matrix_data) |row| {
                allocator.free(row);
            }
            allocator.free(matrix_data);
        }

        for (0..n) |i| {
            matrix_data[i] = try allocator.alloc(f32, n + 1);
        }

        // Initialize the matrix with coefficients
        matrix_data[0][0] = 4.0;  matrix_data[0][1] = -1.0; matrix_data[0][2] = 0.0;  matrix_data[0][3] = 0.0;  matrix_data[0][4] = 0.0;  matrix_data[0][5] = 2.0;
        matrix_data[1][0] = -1.0; matrix_data[1][1] = 4.0;  matrix_data[1][2] = -1.0; matrix_data[1][3] = 0.0;  matrix_data[1][4] = 0.0;  matrix_data[1][5] = 1.0;
        matrix_data[2][0] = 0.0;  matrix_data[2][1] = -1.0; matrix_data[2][2] = 4.0;  matrix_data[2][3] = -1.0; matrix_data[2][4] = 0.0;  matrix_data[2][5] = 1.0;
        matrix_data[3][0] = 0.0;  matrix_data[3][1] = 0.0;  matrix_data[3][2] = -1.0; matrix_data[3][3] = 4.0;  matrix_data[3][4] = -1.0; matrix_data[3][5] = 1.0;
        matrix_data[4][0] = 0.0;  matrix_data[4][1] = 0.0;  matrix_data[4][2] = 0.0;  matrix_data[4][3] = -1.0; matrix_data[4][4] = 4.0;  matrix_data[4][5] = 2.0;

        var matrix = try GaussSeidel.Matrix.init(allocator, matrix_data);
        defer matrix.deinit();

        std.debug.print("Input Matrix:\n", .{});
        matrix.print();

        const solution = try matrix.solve(show_iterations);
        defer allocator.free(solution);
    } else if (std.mem.startsWith(u8, method, "trapezoidal") or
               std.mem.startsWith(u8, method, "simpson") or
               std.mem.startsWith(u8, method, "romberg") or
               std.mem.startsWith(u8, method, "gaussian")) {
        // Run numerical integration methods
        const a: f32 = 0.0;
        const b: f32 = 1.0;
        const n: usize = 100;
        const tolerance: f32 = 0.0001;
        const max_iter: usize = 20;

        if (std.mem.eql(u8, method, "trapezoidal")) {
            const result = try NumericalInteg.trapezoidal(integ1, a, b, n);
            std.debug.print("Trapezoidal Rule result: {d:.6}\n", .{result});
        } else if (std.mem.eql(u8, method, "simpson")) {
            const result = try NumericalInteg.simpson(integ1, a, b, n);
            std.debug.print("Simpson's Rule result: {d:.6}\n", .{result});
        } else if (std.mem.eql(u8, method, "romberg")) {
            const result = try NumericalInteg.romberg(integ1, a, b, tolerance, max_iter);
            std.debug.print("Romberg Integration result: {d:.6}\n", .{result});
        } else if (std.mem.eql(u8, method, "gaussian")) {
            const result = try NumericalInteg.gaussianQuadrature(integ1, a, b);
            std.debug.print("Gaussian Quadrature result: {d:.6}\n", .{result});
        }
    } else if (std.mem.eql(u8, method, "lpp")) {
        // Example LPP:
        // Maximize: 3x1 + 2x2
        // Subject to:
        // 2x1 + x2 <= 100
        // x1 + x2 <= 80
        // x1 <= 40
        // x1, x2 >= 0

        var problem = try LPP.Problem.init(allocator, 2, 3);
        errdefer problem.deinit();

        // Set objective function coefficients
        problem.c[0] = 3.0;
        problem.c[1] = 2.0;

        // Set constraint matrix
        problem.A[0][0] = 2.0; problem.A[0][1] = 1.0;
        problem.A[1][0] = 1.0; problem.A[1][1] = 1.0;
        problem.A[2][0] = 1.0; problem.A[2][1] = 0.0;

        // Set right-hand side values
        problem.b[0] = 100.0;
        problem.b[1] = 80.0;
        problem.b[2] = 40.0;

        // Set constraint senses
        problem.sense[0] = '<';
        problem.sense[1] = '<';
        problem.sense[2] = '<';

        // Set variable types (all continuous)
        problem.var_type[0] = 'C';
        problem.var_type[1] = 'C';

        const solution = try LPP.solve(&problem, show_iterations);
        defer solution.deinit();

        if (!show_iterations) {
            std.debug.print("\nSolution:\n", .{});
            solution.print();
        }
    } else {
        // Run root finding methods
        const tolerance: f32 = 0.0001;
        const max_iter: usize = 100;

        if (std.mem.eql(u8, method, "bisection")) {
            const root = try RootFinding.bisection(f1, 0.0, 3.0, tolerance, max_iter);
            std.debug.print("Bisection root: {d:.6}\n", .{root});
        } else if (std.mem.eql(u8, method, "falsepos")) {
            const root = try RootFinding.falsePosition(f1, 0.0, 3.0, tolerance, max_iter);
            std.debug.print("False Position root: {d:.6}\n", .{root});
        } else if (std.mem.eql(u8, method, "fixedpoint")) {
            const root = try RootFinding.fixedPoint(g1, 1.0, tolerance, max_iter);
            std.debug.print("Fixed Point root: {d:.6}\n", .{root});
        } else if (std.mem.eql(u8, method, "newton")) {
            const root = try RootFinding.newtonRaphson(f2, df2, 1.0, tolerance, max_iter);
            std.debug.print("Newton-Raphson root: {d:.6}\n", .{root});
        } else if (std.mem.eql(u8, method, "brent")) {
            const root = try RootFinding.brent(f2, 1.0, 2.0, tolerance, max_iter);
            std.debug.print("Brent's method root: {d:.6}\n", .{root});
        } else {
            std.debug.print("Unknown method: {s}\n", .{method});
            return;
        }
    }
}
