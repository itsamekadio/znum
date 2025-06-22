const std = @import("std");

/// Root finding algorithms for single-variable functions.
pub const rootfinding = struct {
    /// Errors that can be returned by root finding methods.
    pub const RootFindingError = error{
        /// The function does not change sign in the interval [a, b].
        NoRootInInterval,
        /// The method did not converge within the maximum number of iterations.
        MaxIterationsReached,
        /// The function or its derivative is invalid at some point (e.g., zero derivative in Newton-Raphson).
        InvalidFunction,
    };

    /// Returns the absolute value of a float.
    fn abs(x: f32) f32 {
        return if (x < 0) -x else x;
    }

    /// Finds a root of the function `F` in the interval [a, b] using the bisection method.
    ///
    /// Parameters:
    ///   - F: The function whose root is to be found. Must be continuous on [a, b].
    ///   - a, b: The interval endpoints. Must satisfy F(a) * F(b) < 0.
    ///   - tolerance: The stopping criterion for the root approximation.
    ///   - max_iter: The maximum number of iterations to perform.
    ///
    /// Returns:
    ///   - The approximate root as `f32` on success.
    ///   - `RootFindingError.NoRootInInterval` if F(a) * F(b) >= 0.
    ///   - `RootFindingError.MaxIterationsReached` if the method does not converge within `max_iter`.
    pub fn bisection(comptime F: fn(f32) f32, a: f32, b: f32, tolerance: f32, max_iter: usize) !f32 {
        if (F(a) * F(b) >= 0) return RootFindingError.NoRootInInterval;

        var a_val = a;
        var b_val = b;
        var c: f32 = 0.0;

        for (0..max_iter) |_| {
            c = (a_val + b_val) / 2.0;
            if (F(c) == 0.0 or (b_val - a_val) / 2.0 < tolerance) {
                return c;
            }

            if (F(c) * F(a_val) < 0) {
                b_val = c;
            } else {
                a_val = c;
            }
        }

        return RootFindingError.MaxIterationsReached;
    }

    /// Finds a root of the function `F` in the interval [a, b] using the false position (regula falsi) method.
    ///
    /// Parameters:
    ///   - F: The function whose root is to be found. Must be continuous on [a, b].
    ///   - a, b: The interval endpoints. Must satisfy F(a) * F(b) < 0.
    ///   - tolerance: The stopping criterion for the root approximation.
    ///   - max_iter: The maximum number of iterations to perform.
    ///
    /// Returns:
    ///   - The approximate root as `f32` on success.
    ///   - `RootFindingError.NoRootInInterval` if F(a) * F(b) >= 0.
    ///   - `RootFindingError.MaxIterationsReached` if the method does not converge within `max_iter`.
    pub fn falsePosition(comptime F: fn(f32) f32, a: f32, b: f32, tolerance: f32, max_iter: usize) !f32 {
        if (F(a) * F(b) >= 0) return RootFindingError.NoRootInInterval;

        var a_val = a;
        var b_val = b;
        var c: f32 = 0.0;

        for (0..max_iter) |_| {
            c = (a_val * F(b_val) - b_val * F(a_val)) / (F(b_val) - F(a_val));
            if (F(c) == 0.0 or (b_val - a_val) < tolerance) {
                return c;
            }

            if (F(c) * F(a_val) < 0) {
                b_val = c;
            } else {
                a_val = c;
            }
        }

        return RootFindingError.MaxIterationsReached;
    }

    /// Finds a fixed point of the function `G` using the fixed point iteration method.
    ///
    /// Parameters:
    ///   - G: The iteration function. The root is a value x such that G(x) = x.
    ///   - x0: The initial guess for the root.
    ///   - tolerance: The stopping criterion for the root approximation.
    ///   - max_iter: The maximum number of iterations to perform.
    ///
    /// Returns:
    ///   - The approximate fixed point as `f32` on success.
    ///   - `RootFindingError.MaxIterationsReached` if the method does not converge within `max_iter`.
    pub fn fixedPoint(comptime G: fn(f32) f32, x0: f32, tolerance: f32, max_iter: usize) !f32 {
        var x = x0;
        var x_prev: f32 = 0.0;

        for (0..max_iter) |_| {
            x_prev = x;
            x = G(x);

            if (abs(x - x_prev) < tolerance) {
                return x;
            }
        }

        return RootFindingError.MaxIterationsReached;
    }

    /// Finds a root of the function `F` using the Newton-Raphson method.
    ///
    /// Parameters:
    ///   - F: The function whose root is to be found.
    ///   - DF: The derivative of F.
    ///   - x0: The initial guess for the root.
    ///   - tolerance: The stopping criterion for the root approximation.
    ///   - max_iter: The maximum number of iterations to perform.
    ///
    /// Returns:
    ///   - The approximate root as `f32` on success.
    ///   - `RootFindingError.InvalidFunction` if the derivative is zero at any iteration.
    ///   - `RootFindingError.MaxIterationsReached` if the method does not converge within `max_iter`.
    pub fn newtonRaphson(comptime F: fn(f32) f32, comptime DF: fn(f32) f32, x0: f32, tolerance: f32, max_iter: usize) !f32 {
        var x = x0;

        for (0..max_iter) |_| {
            const fx = F(x);
            const dfx = DF(x);
            
            if (dfx == 0.0) return RootFindingError.InvalidFunction;
            
            const x_next = x - fx / dfx;
            if (abs(x_next - x) < tolerance) {
                return x_next;
            }
            x = x_next;
        }

        return RootFindingError.MaxIterationsReached;
    }

    /// Finds a root of the function `F` in the interval [a, b] using Brent's method.
    ///
    /// Parameters:
    ///   - F: The function whose root is to be found. Must be continuous on [a, b].
    ///   - a, b: The interval endpoints. Must satisfy F(a) * F(b) < 0.
    ///   - tolerance: The stopping criterion for the root approximation.
    ///   - max_iter: The maximum number of iterations to perform.
    ///
    /// Returns:
    ///   - The approximate root as `f32` on success.
    ///   - `RootFindingError.NoRootInInterval` if F(a) * F(b) >= 0.
    ///   - `RootFindingError.MaxIterationsReached` if the method does not converge within `max_iter`.
    pub fn brent(comptime F: fn(f32) f32, a: f32, b: f32, tolerance: f32, max_iter: usize) !f32 {
        if (F(a) * F(b) >= 0) return RootFindingError.NoRootInInterval;

        var a_val = a;
        var b_val = b;
        var c = b_val;
        var d: f32 = 0.0;
        var e: f32 = 0.0;
        var fa = F(a_val);
        var fb = F(b_val);
        var fc = fb;

        for (0..max_iter) |_| {
            if (fb * fc > 0) {
                c = a_val;
                fc = fa;
                d = b_val - a_val;
                e = d;
            }

            if (abs(fc) < abs(fb)) {
                a_val = b_val;
                b_val = c;
                c = a_val;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            const tol1 = 2.0 * std.math.floatEps(f32) * abs(b_val) + 0.5 * tolerance;
            const xm = 0.5 * (c - b_val);

            if (abs(xm) <= tol1 or fb == 0.0) {
                return b_val;
            }

            if (abs(e) >= tol1 and abs(fa) > abs(fb)) {
                var s = fb / fa;
                var p: f32 = 0.0;
                var q: f32 = 0.0;

                if (a_val == c) {
                    // Linear interpolation
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    // Inverse quadratic interpolation
                    const q1 = fa / fc;
                    const r = fb / fc;
                    p = s * (2.0 * xm * q1 * (q1 - r) - (b_val - a_val) * (r - 1.0));
                    q = (q1 - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (p > 0.0) {
                    q = -q;
                } else {
                    p = -p;
                }

                s = e;
                e = d;

                if (2.0 * p < 3.0 * xm * q - abs(tol1 * q) and p < abs(0.5 * s * q)) {
                    d = p / q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;
                e = d;
            }

            a_val = b_val;
            fa = fb;

            if (abs(d) > tol1) {
                b_val += d;
            } else {
                b_val += if (xm > 0) tol1 else -tol1;
            }

            fb = F(b_val);
        }

        return RootFindingError.MaxIterationsReached;
    }
};
