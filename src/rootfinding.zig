const std = @import("std");

pub const Matrix = struct {
    data: [2][3]i32,

    pub fn init(data: [2][3]i32) Matrix {
        return Matrix{ .data = data };
    }

    pub fn print(self: Matrix) void {
        const stdout = std.io.getStdOut().writer();
        for (self.data) |row| {
            for (row) |val| {
                stdout.print("{d} ", .{val}) catch {};
            }
            stdout.print("\n", .{}) catch {};
        }
    }
};

pub const rootfinding = struct {
    pub const RootFindingError = error{
        NoRootInInterval,
        MaxIterationsReached,
        InvalidFunction,
    };

    // Helper function to get absolute value
    fn abs(x: f32) f32 {
        return if (x < 0) -x else x;
    }

    pub fn bisection(f: fn(f32) f32, a: f32, b: f32, tolerance: f32, max_iter: usize) !f32 {
        if (f(a) * f(b) >= 0) return RootFindingError.NoRootInInterval;

        var a_val = a;
        var b_val = b;
        var c: f32 = 0.0;

        for (0..max_iter) |_| {
            c = (a_val + b_val) / 2.0;
            if (f(c) == 0.0 or (b_val - a_val) / 2.0 < tolerance) {
                return c;
            }

            if (f(c) * f(a_val) < 0) {
                b_val = c;
            } else {
                a_val = c;
            }
        }

        return RootFindingError.MaxIterationsReached;
    }

    pub fn falsePosition(f: fn(f32) f32, a: f32, b: f32, tolerance: f32, max_iter: usize) !f32 {
        if (f(a) * f(b) >= 0) return RootFindingError.NoRootInInterval;

        var a_val = a;
        var b_val = b;
        var c: f32 = 0.0;

        for (0..max_iter) |_| {
            c = (a_val * f(b_val) - b_val * f(a_val)) / (f(b_val) - f(a_val));
            if (f(c) == 0.0 or (b_val - a_val) < tolerance) {
                return c;
            }

            if (f(c) * f(a_val) < 0) {
                b_val = c;
            } else {
                a_val = c;
            }
        }

        return RootFindingError.MaxIterationsReached;
    }

    pub fn fixedPoint(g: fn(f32) f32, x0: f32, tolerance: f32, max_iter: usize) !f32 {
        var x = x0;
        var x_prev: f32 = 0.0;

        for (0..max_iter) |_| {
            x_prev = x;
            x = g(x);

            if (abs(x - x_prev) < tolerance) {
                return x;
            }
        }

        return RootFindingError.MaxIterationsReached;
    }

    pub fn newtonRaphson(f: fn(f32) f32, df: fn(f32) f32, x0: f32, tolerance: f32, max_iter: usize) !f32 {
        var x = x0;

        for (0..max_iter) |_| {
            const fx = f(x);
            const dfx = df(x);
            
            if (dfx == 0.0) return RootFindingError.InvalidFunction;
            
            const x_next = x - fx / dfx;
            if (abs(x_next - x) < tolerance) {
                return x_next;
            }
            x = x_next;
        }

        return RootFindingError.MaxIterationsReached;
    }

    pub fn brent(f: fn(f32) f32, a: f32, b: f32, tolerance: f32, max_iter: usize) !f32 {
        if (f(a) * f(b) >= 0) return RootFindingError.NoRootInInterval;

        var a_val = a;
        var b_val = b;
        var c = b_val;
        var d: f32 = 0.0;
        var e: f32 = 0.0;
        var fa = f(a_val);
        var fb = f(b_val);
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

            fb = f(b_val);
        }

        return RootFindingError.MaxIterationsReached;
    }
};
