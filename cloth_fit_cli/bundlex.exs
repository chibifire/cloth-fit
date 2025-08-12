defmodule ClothFitCli.BundlexProject do
  use Bundlex.Project

  def project do
    [
      natives: [
        polyfem_nif: [
          sources: ["polyfem_nif.cpp"],
          interface: :nif,
          preprocessor: Unifex,
          includes: ["../src"],
          libs: ["polyfem"],
          lib_dirs: ["../build"],
          language: :cpp,
          compiler_flags: ["-std=c++17", "-O2"],
          linker_flags: ["-lstdc++"]
        ]
      ]
    ]
  end
end
