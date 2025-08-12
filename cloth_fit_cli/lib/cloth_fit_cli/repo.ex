defmodule ClothFitCli.Repo do
  use Ecto.Repo,
    otp_app: :cloth_fit_cli,
    adapter: Ecto.Adapters.SQLite3
end
