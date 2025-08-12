defmodule ClothFitCli.MixProject do
  use Mix.Project

  def project do
    [
      app: :cloth_fit_cli,
      version: "0.1.0",
      elixir: "~> 1.18",
      start_permanent: Mix.env() == :prod,
      deps: deps()
    ]
  end

  # Run "mix help compile.app" to learn about applications.
  def application do
    [
      extra_applications: [:logger],
      mod: {ClothFitCli.Application, []}
    ]
  end

  # Run "mix help deps" to learn about dependencies.
  defp deps do
    [
      {:oban, "~> 2.15"},
      {:ecto_sqlite3, "~> 0.8"},
      {:phoenix_pubsub, "~> 2.0"},
      {:igniter, "~> 0.6", only: [:dev, :test]}
    ]
  end
end
