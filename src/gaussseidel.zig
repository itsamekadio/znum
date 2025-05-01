const std = @import("std");

pub const gaussseidel = struct {
    pub const Matrix = struct {
        data: [][]f32,
        allocator: std.mem.Allocator,

        pub fn init(allocator: std.mem.Allocator, data: [][]f32) !Matrix {
            return Matrix{
                .data = data,
                .allocator = allocator,
            };
        }

        pub fn deinit(self: Matrix) void {
            for (self.data) |row| {
                self.allocator.free(row);
            }
            self.allocator.free(self.data);
        }

        pub fn print(self: Matrix) void {
            const stdout = std.io.getStdOut().writer();
            for (self.data) |row| {
                for (row) |val| {
                    stdout.print("{d:.2} ", .{val}) catch {};
                }
                stdout.print("\n", .{}) catch {};
            }
        }

        pub fn solve(self: Matrix, show_iterations: bool) ![]f32 {
            const n = self.data.len;
            const max_iter = 25;
            const tolerance = 0.001;

            var x = try self.allocator.alloc(f32, n);
            var prev = try self.allocator.alloc(f32, n);
            defer self.allocator.free(prev);
            
            @memset(x, 0.0);
            @memset(prev, 0.0);

            const stdout = std.io.getStdOut().writer();

            for (0..max_iter) |iter| {
                for (0..n) |i| {
                    var sum: f32 = self.data[i][n];
                    for (0..n) |j| {
                        if (i != j) {
                            sum -= self.data[i][j] * x[j];
                        }
                    }
                    x[i] = sum / self.data[i][i];
                }

                if (show_iterations) {
                    stdout.print("Iteration {d}: ", .{iter + 1}) catch {};
                    for (x) |xi| {
                        stdout.print("{d:.4} ", .{xi}) catch {};
                    }
                    stdout.print("\n", .{}) catch {};
                }

                // Convergence check
                var converged = true;
                for (0..n) |i| {
                    var diff = x[i] - prev[i];
                    if (diff < 0) {
                        diff = -diff;
                    }
                    
                    if (diff > tolerance) {
                        converged = false;
                        break;
                    }
                }

                if (converged) {
                    stdout.print("Converged at iteration {d}\n", .{iter + 1}) catch {};
                    break;
                }

                // Update prev with current x values
                for (0..n) |i| {
                    prev[i] = x[i];
                }
            }

            stdout.print("Final Solution:\n", .{}) catch {};
            for (x) |xi| {
                stdout.print("{d:.4}\n", .{xi}) catch {};
            }

            return x;
        }
    };
};
