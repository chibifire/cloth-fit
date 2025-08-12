defmodule ClothFitCli.MixProject do
  use Mix.Project

  def project do
    [
      app: :cloth_fit_cli,
      version: "0.1.0",
      elixir: "~> 1.18",
      start_permanent: Mix.env() == :prod,
      deps: deps(),
      compilers: [:bundlex] ++ Mix.compilers(),
      aliases: aliases()
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
      {:unifex, "~> 1.0"},
      {:bundlex, "~> 1.0"},
      {:ecto, "~> 3.12"},
      {:ecto_sql, "~> 3.12"},
      {:ecto_sqlite3, "~> 0.17"},
      {:oban, "~> 2.18"},
      {:jason, "~> 1.4"},
      {:igniter, "~> 0.6", only: [:dev, :test]}
    ]
  end

  defp aliases do
    [
      setup: ["deps.get", "ecto.setup"],
      "ecto.setup": ["ecto.create", "ecto.migrate", "run priv/repo/seeds.exs"],
      "ecto.reset": ["ecto.drop", "ecto.setup"],
      test: ["ecto.create --quiet", "ecto.migrate --quiet", "test"]
    ]
  end
end
