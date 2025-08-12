module PolyFem

interface [NIF, CNode]

spec simulate(config :: payload, output_path :: string) :: {:ok :: label, result :: string} | {:error :: label, reason :: string}
spec validate_garment_mesh(mesh_path :: string) :: {:ok :: label, valid :: bool} | {:error :: label, reason :: string}
spec validate_avatar_mesh(mesh_path :: string) :: {:ok :: label, valid :: bool} | {:error :: label, reason :: string}
spec load_garment_info(garment_path :: string) :: {:ok :: label, info :: string} | {:error :: label, reason :: string}
spec load_avatar_info(avatar_path :: string) :: {:ok :: label, info :: string} | {:error :: label, reason :: string}
