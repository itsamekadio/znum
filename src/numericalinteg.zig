const std = @import("std");

/// Numerical integration algorithms for single-variable functions.
pub const numericalinteg = struct {
    /// Errors that can be returned by integration methods.
    pub const IntegrationError = error{
        /// The interval [a, b] is invalid (a >= b).
        InvalidInterval,
        /// The function or parameters are invalid (e.g., n = 0, n not divisible by required value).
        InvalidFunction,
        /// The method did not converge within the maximum number of iterations (not used here).
        MaxIterationsReached,
    };

    /// Returns the absolute value of a float.
    fn abs(x: f32) f32 {
        return if (x < 0) -x else x;
    }

    /// Approximates the definite integral of `F` over [a, b] using the trapezoidal rule.
    ///
    /// Parameters:
    ///   - F: The function to integrate.
    ///   - a, b: The interval endpoints (a < b).
    ///   - n: The number of subintervals (n > 0).
    ///
    /// Returns:
    ///   - The approximate integral as `f32` on success.
    ///   - `IntegrationError.InvalidInterval` if a >= b.
    ///   - `IntegrationError.InvalidFunction` if n == 0.
    pub fn trapezoidal(comptime F: fn(f32) f32, a: f32, b: f32, n: usize) !f32 {
        if (a >= b) return IntegrationError.InvalidInterval;
        if (n == 0) return IntegrationError.InvalidFunction;

        const h = (b - a) / @as(f32, @floatFromInt(n));
        var sum: f32 = 0.5 * (F(a) + F(b));

        var i: usize = 1;
        while (i < n) : (i += 1) {
            const x = a + @as(f32, @floatFromInt(i)) * h;
            sum += F(x);
        }

        return sum * h;
    }

    /// Approximates the definite integral of `F` over [a, b] using Simpson's 1/3 rule.
    ///
    /// Parameters:
    ///   - F: The function to integrate.
    ///   - a, b: The interval endpoints (a < b).
    ///   - n: The number of subintervals (must be even and > 0).
    ///
    /// Returns:
    ///   - The approximate integral as `f32` on success.
    ///   - `IntegrationError.InvalidInterval` if a >= b.
    ///   - `IntegrationError.InvalidFunction` if n == 0 or n is not even.
    pub fn simpson13(comptime F: fn(f32) f32, a: f32, b: f32, n: usize) !f32 {
        if (a >= b) return IntegrationError.InvalidInterval;
        if (n == 0 or n % 2 != 0) return IntegrationError.InvalidFunction;

        const h = (b - a) / @as(f32, @floatFromInt(n));
        var sum: f32 = F(a) + F(b);

        var i: usize = 1;
        while (i < n) : (i += 1) {
            const x = a + @as(f32, @floatFromInt(i)) * h;
            sum += if (i % 2 == 0) 2.0 * F(x) else 4.0 * F(x);
        }

        return sum * h / 3.0;
    }

    /// Approximates the definite integral of `F` over [a, b] using Simpson's 3/8 rule.
    ///
    /// Parameters:
    ///   - F: The function to integrate.
    ///   - a, b: The interval endpoints (a < b).
    ///   - n: The number of subintervals (must be a multiple of 3 and > 0).
    ///
    /// Returns:
    ///   - The approximate integral as `f32` on success.
    ///   - `IntegrationError.InvalidInterval` if a >= b.
    ///   - `IntegrationError.InvalidFunction` if n == 0 or n is not a multiple of 3.
    pub fn simpson38(comptime F: fn(f32) f32, a: f32, b: f32, n: usize) !f32 {
        if (a >= b) return IntegrationError.InvalidInterval;
        if (n == 0 or n % 3 != 0) return IntegrationError.InvalidFunction;

        const h = (b - a) / @as(f32, @floatFromInt(n));
        var sum: f32 = F(a) + F(b);

        var i: usize = 1;
        while (i < n) : (i += 1) {
            const x = a + @as(f32, @floatFromInt(i)) * h;
            sum += if (i % 3 == 0) 2.0 * F(x) else 3.0 * F(x);
        }

        return sum * 3.0 * h / 8.0;
    }

    /// Approximates the definite integral of `F` over [a, b] using 2-point Gaussian quadrature.
    ///
    /// Parameters:
    ///   - F: The function to integrate.
    ///   - a, b: The interval endpoints (a < b).
    ///
    /// Returns:
    ///   - The approximate integral as `f32` on success.
    ///   - `IntegrationError.InvalidInterval` if a >= b.
    pub fn gaussianQuadrature(comptime F: fn(f32) f32, a: f32, b: f32) !f32 {
        if (a >= b) return IntegrationError.InvalidInterval;

        // 2-point Gaussian quadrature weights and points (stack allocated)
        const weights = [_]f32{ 1.0, 1.0 };
        const points = [_]f32{ -0.5773502691896257, 0.5773502691896257 };

        // Transform points from [-1,1] to [a,b]
        const c1 = (b - a) / 2.0;
        const c2 = (b + a) / 2.0;

        var sum: f32 = 0.0;
        for (points, 0..) |x, i| {
            const transformed_x = c1 * x + c2;
            sum += weights[i] * F(transformed_x);
        }

        return c1 * sum;
    }
};
