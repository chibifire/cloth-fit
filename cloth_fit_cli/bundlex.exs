defmodule PolyFem.BundlexProject do
  use Bundlex.Project

  def project() do
  [
    natives: natives(Bundlex.get_target())
  ]
  end

  def natives(_platform) do
  [
      polyfem: [
        sources: ["polyfem.cpp"],
        interface: [:nif, :cnode],
        preprocessor: Unifex,
        language: :cpp,
        includes: ["../../src"],
        compiler_flags: ["-std=c++17", "-O2"],
        linker_flags: ["-lstdc++"]
      ]
  ]
  end
end
