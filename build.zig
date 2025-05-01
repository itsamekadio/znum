const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Define the matrix module
    const matrix_mod = b.createModule(.{
        .root_source_file = b.path("src/rootfinding.zig"),
        .target = target,
        .optimize = optimize,
    });
    const numericalinteg_mod = b.createModule(.{
        .root_source_file = b.path("src/numericalinteg.zig"),
        .target = target,
        .optimize = optimize,
    });
    const lpp_mod = b.createModule(.{
        .root_source_file = b.path("src/lpp.zig"),
        .target = target,
        .optimize = optimize,
    });
     const simplex_mod = b.createModule(.{
        .root_source_file = b.path("src/gaussseidel.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Define the library module
    const lib_mod = b.createModule(.{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Define the executable module
    const exe_mod = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Add imports
    exe_mod.addImport("proj_lib", lib_mod);
    exe_mod.addImport("matrix", matrix_mod);
    exe_mod.addImport("simplex", simplex_mod);
    exe_mod.addImport("numericalinteg", numericalinteg_mod);
    exe_mod.addImport("lpp", lpp_mod);




    // Create the static library
    const lib = b.addLibrary(.{
        .linkage = .static,
        .name = "proj",
        .root_module = lib_mod,
    });

    // Install the library
    b.installArtifact(lib);

    // Create the executable
    const exe = b.addExecutable(.{
        .name = "proj",
        .root_module = exe_mod,
    });

    // Install the executable
    b.installArtifact(exe);

    // Add run command
    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());

    // Pass command-line arguments to the run command
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    // Define the run step
    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    // Define unit tests for the library
    const lib_unit_tests = b.addTest(.{
        .root_module = lib_mod,
    });

    // Define the run command for library unit tests
    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);

    // Define unit tests for the executable
    const exe_unit_tests = b.addTest(.{
        .root_module = exe_mod,
    });

    // Define the run command for executable unit tests
    const run_exe_unit_tests = b.addRunArtifact(exe_unit_tests);

    // Define the test step
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);
    test_step.dependOn(&run_exe_unit_tests.step);
}
