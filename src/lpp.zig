const std = @import("std");

pub const lpp = struct {
    pub const LPPError = error{
        InvalidProblem,
        Infeasible,
        Unbounded,
        MaxIterationsReached,
    };

    pub const Problem = struct {
        allocator: std.mem.Allocator,
        c: []f32, // Objective function coefficients
        A: [][]f32, // Constraint matrix
        b: []f32, // Right-hand side values
        sense: []u8, // Constraint senses ("<=", ">=", "=")
        var_type: []u8, // Variable types ("C" for continuous, "I" for integer)
        n_vars: usize,
        n_constraints: usize,

        pub fn init(allocator: std.mem.Allocator, n_vars: usize, n_constraints: usize) !Problem {
            const c = try allocator.alloc(f32, n_vars);
            errdefer allocator.free(c);

            const A = try allocator.alloc([]f32, n_constraints);
            errdefer {
                for (A) |row| {
                    allocator.free(row);
                }
                allocator.free(A);
            }

            for (0..n_constraints) |i| {
                A[i] = try allocator.alloc(f32, n_vars);
            }

            const b = try allocator.alloc(f32, n_constraints);
            errdefer allocator.free(b);

            const sense = try allocator.alloc(u8, n_constraints);
            errdefer allocator.free(sense);

            const var_type = try allocator.alloc(u8, n_vars);
            errdefer allocator.free(var_type);

            return Problem{
                .allocator = allocator,
                .c = c,
                .A = A,
                .b = b,
                .sense = sense,
                .var_type = var_type,
                .n_vars = n_vars,
                .n_constraints = n_constraints,
            };
        }

        pub fn deinit(self: *const Problem) void {
            for (self.A) |row| {
                self.allocator.free(row);
            }
            self.allocator.free(self.A);
            self.allocator.free(self.c);
            self.allocator.free(self.b);
            self.allocator.free(self.sense);
            self.allocator.free(self.var_type);
        }

        pub fn print(self: *const Problem) void {
            std.debug.print("Linear Programming Problem:\n", .{});
            std.debug.print("Objective: max ", .{});
            for (0..self.n_vars) |i| {
                if (i > 0 and self.c[i] >= 0) std.debug.print("+", .{});
                std.debug.print("{d}x{} ", .{ self.c[i], i + 1 });
            }
            std.debug.print("\n", .{});

            std.debug.print("Subject to:\n", .{});
            for (0..self.n_constraints) |i| {
                for (0..self.n_vars) |j| {
                    if (j > 0 and self.A[i][j] >= 0) std.debug.print("+", .{});
                    std.debug.print("{d}x{} ", .{ self.A[i][j], j + 1 });
                }
                std.debug.print("{c} {d}\n", .{ self.sense[i], self.b[i] });
            }
        }
    };

    pub const Solution = struct {
        allocator: std.mem.Allocator,
        x: []f32, // Variable values
        obj_value: f32, // Objective value
        status: []const u8, // Solution status

        pub fn init(allocator: std.mem.Allocator, n_vars: usize) !Solution {
            const x = try allocator.alloc(f32, n_vars);
            return Solution{
                .allocator = allocator,
                .x = x,
                .obj_value = 0.0,
                .status = "Optimal",
            };
        }

        pub fn deinit(self: *const Solution) void {
            self.allocator.free(self.x);
        }

        pub fn print(self: *const Solution) void {
            std.debug.print("Solution status: {s}\n", .{self.status});
            std.debug.print("Objective value: {d}\n", .{self.obj_value});
            std.debug.print("Variable values:\n", .{});
            for (self.x, 0..) |val, i| {
                std.debug.print("x{} = {d}\n", .{ i + 1, val });
            }
        }
    };

    pub fn solve(problem: *Problem, show_iterations: bool) !Solution {
        if (show_iterations) {
            std.debug.print("\nInitial Problem Formulation:\n", .{});
            problem.print();
        }

        // Convert to standard form
        const std_form = try convertToStandardForm(problem);
        defer std_form.deinit();

        if (show_iterations) {
            std.debug.print("\nStandard Form:\n", .{});
            std_form.print();
        }

        // Initialize tableau
        const tableau = try createTableau(problem);
        defer {
            for (tableau) |row| {
                problem.allocator.free(row);
            }
            problem.allocator.free(tableau);
        }

        if (show_iterations) {
            std.debug.print("\nInitial Tableau:\n", .{});
            printTableau(tableau);
        }

        // Perform simplex iterations
        var solution = try Solution.init(problem.allocator, problem.n_vars);
        errdefer solution.deinit();

        var iter: usize = 0;
        const max_iter: usize = 1000;

        while (iter < max_iter) : (iter += 1) {
            // Find entering variable (most negative reduced cost)
            const entering = findEnteringVariable(tableau);
            if (entering == -1) {
                // Optimal solution found
                break;
            }

            // Find leaving variable (minimum ratio test)
            const leaving = findLeavingVariable(tableau, entering);
            if (leaving == -1) {
                solution.status = "Unbounded";
                return solution;
            }

            if (show_iterations) {
                std.debug.print("\nIteration {}:\n", .{iter + 1});
                std.debug.print("Entering variable: x{}\n", .{entering + 1});
                std.debug.print("Leaving variable: x{}\n", .{leaving + 1});
            }

            // Pivot
            performPivot(tableau, entering, leaving);

            if (show_iterations) {
                std.debug.print("Tableau after pivot:\n", .{});
                printTableau(tableau);
            }
        }

        if (iter >= max_iter) {
            return LPPError.MaxIterationsReached;
        }

        // Extract solution
        try extractSolution(tableau, &solution);

        return solution;
    }

    fn printTableau(tableau: [][]f32) void {
        const n_rows = tableau.len;
        const n_cols = tableau[0].len;

        // Print column headers
        std.debug.print("Basis | ", .{});
        for (0..n_cols-1) |j| {
            std.debug.print("x{} | ", .{j + 1});
        }
        std.debug.print("RHS\n", .{});

        // Print separator
        for (0..n_cols) |_| {
            std.debug.print("------", .{});
        }
        std.debug.print("\n", .{});

        // Print rows
        for (0..n_rows) |i| {
            if (i == 0) {
                std.debug.print("Obj  | ", .{});
            } else {
                std.debug.print("x{}   | ", .{i});
            }
            for (0..n_cols) |j| {
                std.debug.print("{d:6.2} | ", .{tableau[i][j]});
            }
            std.debug.print("\n", .{});
        }
    }

    fn convertToStandardForm(problem: *Problem) !Problem {
        // TODO: Implement conversion to standard form
        // This is a placeholder that just returns the original problem
        return problem.*;
    }

    fn createTableau(problem: *const Problem) ![][]f32 {
        const n_rows = problem.n_constraints + 1;
        const n_cols = problem.n_vars + problem.n_constraints + 1;
        
        var tableau = try problem.allocator.alloc([]f32, n_rows);
        errdefer {
            for (tableau) |row| {
                problem.allocator.free(row);
            }
            problem.allocator.free(tableau);
        }

        for (0..n_rows) |i| {
            tableau[i] = try problem.allocator.alloc(f32, n_cols);
        }

        // Initialize tableau with problem data
        // First row: objective function
        for (0..problem.n_vars) |j| {
            tableau[0][j] = -problem.c[j]; // Negative because we're maximizing
        }
        for (problem.n_vars..n_cols-1) |j| {
            tableau[0][j] = 0.0;
        }
        tableau[0][n_cols-1] = 0.0; // RHS

        // Constraint rows
        for (0..problem.n_constraints) |i| {
            for (0..problem.n_vars) |j| {
                tableau[i+1][j] = problem.A[i][j];
            }
            for (problem.n_vars..n_cols-1) |j| {
                tableau[i+1][j] = if (j - problem.n_vars == i) 1.0 else 0.0;
            }
            tableau[i+1][n_cols-1] = problem.b[i];
        }

        return tableau;
    }

    fn findEnteringVariable(tableau: [][]f32) i32 {
        const n_cols = tableau[0].len;
        var min_col: i32 = -1;
        var min_val: f32 = 0.0;

        for (0..n_cols-1) |j| {
            if (tableau[0][j] < min_val) {
                min_val = tableau[0][j];
                min_col = @as(i32, @intCast(j));
            }
        }

        return min_col;
    }

    fn findLeavingVariable(tableau: [][]f32, entering: i32) i32 {
        const n_rows = tableau.len;
        const n_cols = tableau[0].len;
        var min_row: i32 = -1;
        var min_ratio: f32 = std.math.inf(f32);

        for (1..n_rows) |i| {
            if (tableau[i][@intCast(entering)] > 0) {
                const ratio = tableau[i][n_cols-1] / tableau[i][@intCast(entering)];
                if (ratio < min_ratio) {
                    min_ratio = ratio;
                    min_row = @as(i32, @intCast(i));
                }
            }
        }

        return min_row;
    }

    fn performPivot(tableau: [][]f32, entering: i32, leaving: i32) void {
        const n_rows = tableau.len;
        const n_cols = tableau[0].len;
        const pivot_row = @as(usize, @intCast(leaving));
        const pivot_col = @as(usize, @intCast(entering));
        const pivot_val = tableau[pivot_row][pivot_col];

        // Normalize pivot row
        for (0..n_cols) |j| {
            tableau[pivot_row][j] /= pivot_val;
        }

        // Update other rows
        for (0..n_rows) |i| {
            if (i != pivot_row) {
                const factor = tableau[i][pivot_col];
                for (0..n_cols) |j| {
                    tableau[i][j] -= factor * tableau[pivot_row][j];
                }
            }
        }
    }

    fn extractSolution(tableau: [][]f32, solution: *Solution) !void {
        const n_rows = tableau.len;
        const n_cols = tableau[0].len;
        const n_vars = solution.x.len;

        // Extract objective value
        solution.obj_value = tableau[0][n_cols-1];

        // Extract variable values
        for (0..n_vars) |j| {
            var found = false;
            for (1..n_rows) |i| {
                if (tableau[i][j] == 1.0) {
                    var is_basic = true;
                    for (0..n_vars) |k| {
                        if (k != j and tableau[i][k] != 0.0) {
                            is_basic = false;
                            break;
                        }
                    }
                    if (is_basic) {
                        solution.x[j] = tableau[i][n_cols-1];
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                solution.x[j] = 0.0;
            }
        }
    }
};
