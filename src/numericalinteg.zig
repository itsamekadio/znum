const std = @import("std");

pub const numericalinteg = struct {
    pub const IntegrationError = error{
        InvalidInterval,
        InvalidFunction,
        MaxIterationsReached,
    };

    // Helper function to get absolute value
    fn abs(x: f32) f32 {
        return if (x < 0) -x else x;
    }

    // Trapezoidal Rule
    pub fn trapezoidal(f: fn(f32) f32, a: f32, b: f32, n: usize) !f32 {
        if (a >= b) return IntegrationError.InvalidInterval;
        if (n == 0) return IntegrationError.InvalidFunction;

        const h = (b - a) / @as(f32, @floatFromInt(n));
        var sum: f32 = 0.5 * (f(a) + f(b));

        var i: usize = 1;
        while (i < n) : (i += 1) {
            const x = a + @as(f32, @floatFromInt(i)) * h;
            sum += f(x);
        }

        return sum * h;
    }

    // Simpson's Rule
    pub fn simpson(f: fn(f32) f32, a: f32, b: f32, n: usize) !f32 {
        if (a >= b) return IntegrationError.InvalidInterval;
        if (n % 2 != 0) return IntegrationError.InvalidFunction;

        const h = (b - a) / @as(f32, @floatFromInt(n));
        var sum: f32 = f(a) + f(b);

        var i: usize = 1;
        while (i < n) : (i += 1) {
            const x = a + @as(f32, @floatFromInt(i)) * h;
            sum += if (i % 2 == 0) 2.0 * f(x) else 4.0 * f(x);
        }

        return sum * h / 3.0;
    }

    // Romberg Integration
    pub fn romberg(f: fn(f32) f32, a: f32, b: f32, tolerance: f32, max_iter: usize) !f32 {
        if (a >= b) return IntegrationError.InvalidInterval;

        var r = try std.heap.page_allocator.alloc([]f32, max_iter);
        defer {
            for (r) |row| {
                std.heap.page_allocator.free(row);
            }
            std.heap.page_allocator.free(r);
        }

        for (0..max_iter) |i| {
            r[i] = try std.heap.page_allocator.alloc(f32, i + 1);
        }

        // First column using trapezoidal rule
        r[0][0] = 0.5 * (b - a) * (f(a) + f(b));

        for (1..max_iter) |i| {
            // Compute trapezoidal rule for this level
            const power_of_two = @as(usize, 1) << @as(u6, @intCast(i));
            const h = (b - a) / @as(f32, @floatFromInt(power_of_two));
            var sum: f32 = 0.0;
            const max_k = @as(usize, 1) << @as(u6, @intCast(i - 1));
            var k: usize = 1;
            while (k <= max_k) : (k += 1) {
                const x = a + (2.0 * @as(f32, @floatFromInt(k)) - 1.0) * h;
                sum += f(x);
            }
            r[i][0] = 0.5 * r[i-1][0] + h * sum;

            // Richardson extrapolation
            for (1..i+1) |j| {
                const pow4j = std.math.pow(f32, 4.0, @as(f32, @floatFromInt(j)));
                r[i][j] = r[i][j-1] + (r[i][j-1] - r[i-1][j-1]) / (pow4j - 1.0);
            }

            // Check convergence
            if (i > 0 and abs(r[i][i] - r[i-1][i-1]) < tolerance) {
                return r[i][i];
            }
        }

        return IntegrationError.MaxIterationsReached;
    }

    // Gaussian Quadrature (2-point)
    pub fn gaussianQuadrature(f: fn(f32) f32, a: f32, b: f32) !f32 {
        if (a >= b) return IntegrationError.InvalidInterval;

        // 2-point Gaussian quadrature weights and points
        const weights = [_]f32{ 1.0, 1.0 };
        const points = [_]f32{ -0.5773502691896257, 0.5773502691896257 };

        // Transform points from [-1,1] to [a,b]
        const c1 = (b - a) / 2.0;
        const c2 = (b + a) / 2.0;

        var sum: f32 = 0.0;
        for (points, 0..) |x, i| {
            const transformed_x = c1 * x + c2;
            sum += weights[i] * f(transformed_x);
        }

        return c1 * sum;
    }
};
