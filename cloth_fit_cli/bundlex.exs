defmodule PolyFem.BundlexProject do
  use Bundlex.Project

  def project() do
  [
    natives: natives(Bundlex.get_target())
  ]
  end

  def natives(_platform) do
    # Get the absolute path to the build directory
    build_dir = Path.expand("../../_build", __DIR__)

    # Core PolyFEM libraries
    polyfem_libs = [
      "#{build_dir}/libpolyfem.a"
    ]

    # PolyFEM dependencies
    polyfem_deps = [
      "#{build_dir}/_deps/polysolve-build/libpolysolve.a",
      "#{build_dir}/_deps/ipc-toolkit-build/libipc_toolkit.a",
      "#{build_dir}/_deps/tight-inclusion-build/libtight_inclusion.a",
      "#{build_dir}/_deps/scalable-ccd-build/libscalable_ccd.a",
      "#{build_dir}/_deps/finite-diff-build/libfinitediff_finitediff.a",
      "#{build_dir}/_deps/json-spec-engine-build/libjse.a"
    ]

    # SuiteSparse libraries
    suitesparse_libs = [
      "#{build_dir}/libcholmod.a",
      "#{build_dir}/libcolamd.a",
      "#{build_dir}/libamd.a",
      "#{build_dir}/libcamd.a",
      "#{build_dir}/libccolamd.a",
      "#{build_dir}/libsuitesparseconfig.a"
    ]

    # Threading and utility libraries
    threading_libs = [
      "#{build_dir}/_deps/onetbb-build/src/tbb/libtbb.a"
    ]

    # Logging and utility libraries
    utility_libs = [
      "#{build_dir}/_deps/spdlog-build/libspdlog.a",
      "#{build_dir}/_deps/filib-build/libfilib_filib.a"
    ]

    # OpenVDB and related libraries
    openvdb_libs = [
      "#{build_dir}/_deps/openvdb-build/openvdb/openvdb/libopenvdb.a"
    ]

    # Abseil libraries (required by some dependencies)
    abseil_libs = [
      "#{build_dir}/_deps/abseil-cpp-build/absl/base/libabsl_base.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/base/libabsl_spinlock_wait.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/base/libabsl_log_severity.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/base/libabsl_raw_logging_internal.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/base/libabsl_throw_delegate.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/strings/libabsl_strings.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/strings/libabsl_strings_internal.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/numeric/libabsl_int128.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/hash/libabsl_hash.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/hash/libabsl_city.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/hash/libabsl_low_level_hash.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/types/libabsl_bad_optional_access.a",
      "#{build_dir}/_deps/abseil-cpp-build/absl/types/libabsl_bad_variant_access.a"
    ]

    # Combine all libraries in proper linking order
    all_libs = polyfem_libs ++ polyfem_deps ++ suitesparse_libs ++
               threading_libs ++ utility_libs ++ openvdb_libs ++ abseil_libs

    [
      polyfem: [
        sources: ["polyfem.cpp"],
        interface: [:nif, :cnode],
        preprocessor: Unifex,
        language: :cpp,
        includes: [
          Path.expand("../src", __DIR__),
          "#{build_dir}/_deps/eigen-src",
          "#{build_dir}/_deps/spdlog-src/include",
          "#{build_dir}/_deps/nlohmann_json-src/include",
          "#{build_dir}/_deps/ipc-toolkit-src/src",
          "#{build_dir}/_deps/polysolve-src/src",
          "#{build_dir}/_deps/openvdb-src/openvdb",
          "#{build_dir}/_deps/libigl-src/include"
        ],
        compiler_flags: [
          "-std=c++17",
          "-O2",
          "-DEIGEN_STACK_ALLOCATION_LIMIT=0",
          "-DPOLYFEM_SMALL_N=80",
          "-DPOLYFEM_BIG_N=1000",
          "-DSPDLOG_COMPILED_LIB"
        ],
        linker_flags: ["-lstdc++", "-lm", "-lpthread", "-ldl"] ++ all_libs
      ]
    ]
  end
end
