# Numerical Methods in Zig

A collection of numerical methods implemented in Zig, providing efficient and type-safe solutions for:
- Linear System Solvers (Gauss-Seidel)
- Root Finding Methods
- Numerical Integration
- Linear Programming (Simplex Method)

## Quick Start

### Installation
```bash
git clone https://github.com/itsamekadio/numerical-methods-zig
cd numerical-methods-zig
zig build
```

### Basic Usage

#### 1. Linear Systems (Gauss-Seidel)
```bash
# Solve a system of equations
zig build run -- --method gaussseidel

# Show iteration progress
zig build run -- --method gaussseidel --show-iterations
```

#### 2. Root Finding
```bash
# Newton-Raphson method
zig build run -- --method newton

# Bisection method
zig build run -- --method bisection

# Brent's method
zig build run -- --method brent
```

#### 3. Numerical Integration
```bash
# Simpson's Rule
zig build run -- --method simpson

# Trapezoidal Rule
zig build run -- --method trapezoidal

# Romberg Integration
zig build run -- --method romberg
```

#### 4. Linear Programming
```bash
# Solve LP problem
zig build run -- --method lpp

# Show simplex iterations
zig build run -- --method lpp --show-iterations
```

## For Developers

### Project Structure
```
src/
├── main.zig          # Entry point and examples
├── gaussseidel.zig   # Linear system solver
├── rootfinding.zig   # Root finding methods
├── numericalinteg.zig # Integration methods
└── lpp.zig          # Linear programming solver
```

### Adding a New Problem

#### Linear Systems
```zig
const std = @import("std");
const GaussSeidel = @import("gaussseidel.zig");

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    // Create 3x3 system: Ax = b
    const n = 3;
    var matrix = try allocator.alloc([]f32, n);
    defer {
        for (matrix) |row| allocator.free(row);
        allocator.free(matrix);
    }

    // Initialize matrix
    for (0..n) |i| {
        matrix[i] = try allocator.alloc(f32, n + 1);
    }

    // Set coefficients (example system)
    matrix[0][0] = 4.0; matrix[0][1] = -1.0; matrix[0][2] = 0.0; matrix[0][3] = 7.0;
    matrix[1][0] = -1.0; matrix[1][1] = 4.0; matrix[1][2] = -1.0; matrix[1][3] = 6.0;
    matrix[2][0] = 0.0; matrix[2][1] = -1.0; matrix[2][2] = 4.0; matrix[2][3] = 5.0;

    // Solve
    const solution = try GaussSeidel.solve(matrix, 1e-6, 100, true);
    defer allocator.free(solution);

    // Print solution
    for (solution, 0..) |x, i| {
        std.debug.print("x{} = {d:.6}\n", .{ i + 1, x });
    }
}
```

#### Root Finding
```zig
const std = @import("std");
const RootFinding = @import("rootfinding.zig");

// Define function and derivative
fn f(x: f32) f32 {
    return x * x - 4.0; // f(x) = x² - 4
}

fn df(x: f32) f32 {
    return 2.0 * x; // f'(x) = 2x
}

pub fn main() !void {
    // Newton-Raphson
    const root = try RootFinding.newtonRaphson(
        f,          // function
        df,         // derivative
        2.0,        // initial guess
        1e-6,       // tolerance
        100         // max iterations
    );
    std.debug.print("Root found: {d:.6}\n", .{root});
}
```

#### Numerical Integration
```zig
const std = @import("std");
const NumericalInteg = @import("numericalinteg.zig");

// Define function to integrate
fn f(x: f32) f32 {
    return x * x; // f(x) = x²
}

pub fn main() !void {
    // Simpson's Rule
    const integral = try NumericalInteg.simpson(
        f,          // function
        0.0,        // lower limit
        1.0,        // upper limit
        100         // number of intervals
    );
    std.debug.print("Integral: {d:.6}\n", .{integral});
}
```

#### Linear Programming
```zig
const std = @import("std");
const LPP = @import("lpp.zig");

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    // Create LP problem
    var problem = try LPP.Problem.init(allocator, 2, 3);
    defer problem.deinit();

    // Set objective: maximize 3x1 + 2x2
    problem.c[0] = 3.0;
    problem.c[1] = 2.0;

    // Set constraints
    problem.A[0][0] = 2.0; problem.A[0][1] = 1.0; problem.b[0] = 100.0;
    problem.A[1][0] = 1.0; problem.A[1][1] = 1.0; problem.b[1] = 80.0;
    problem.A[2][0] = 1.0; problem.A[2][1] = 0.0; problem.b[2] = 40.0;

    // Set constraint senses (≤)
    problem.sense[0] = '<';
    problem.sense[1] = '<';
    problem.sense[2] = '<';

    // Solve
    const solution = try LPP.solve(&problem, true);
    defer solution.deinit();

    // Print solution
    std.debug.print("Optimal solution:\n", .{});
    for (solution.x, 0..) |value, i| {
        std.debug.print("x{} = {d:.2}\n", .{ i + 1, value });
    }
    std.debug.print("Objective value: {d:.2}\n", .{solution.obj_value});
}
```

### Error Handling
```zig
// Example error handling
const result = LPP.solve(&problem, true) catch |err| {
    switch (err) {
        error.Infeasible => std.debug.print("Problem is infeasible\n", .{}),
        error.Unbounded => std.debug.print("Problem is unbounded\n", .{}),
        error.MaxIterations => std.debug.print("Maximum iterations reached\n", .{}),
        else => std.debug.print("Unknown error: {}\n", .{err}),
    }
    return;
};
```

### Memory Management
```zig
// Example memory management
var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
defer arena.deinit();
const allocator = arena.allocator();

// Allocate with cleanup
var matrix = try allocator.alloc([]f32, n);
errdefer {
    for (matrix) |row| allocator.free(row);
    allocator.free(matrix);
}
```

## Testing

### Run Tests
```bash
zig test src/main.zig
zig test src/gaussseidel.zig
zig test src/rootfinding.zig
zig test src/numericalinteg.zig
zig test src/lpp.zig
```

### Benchmarking
```bash
# Build in release mode
zig build -Doptimize=ReleaseFast

# Run benchmarks
zig build run -- --method gaussseidel --benchmark
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## License

MIT License - see LICENSE file for details