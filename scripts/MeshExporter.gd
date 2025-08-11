class_name MeshExporter
extends RefCounted

const Utils = preload("Utils.gd")

func export_meshes_as_obj(node: Node, output_path: String, avatar_transform: Dictionary = {}):
    print("DEBUG: MeshExporter received transform: ", avatar_transform)
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if not file:
        print("Error creating output file: ", output_path)
        return
    
    file.store_line("# OBJ file exported from Godot")
    file.store_line("mtllib " + output_path.get_basename() + ".mtl")
    var vertex_count = 1
    
    vertex_count = export_node_meshes(node, file, vertex_count, avatar_transform)
    file.close()

func export_node_meshes(node: Node, file: FileAccess, vertex_count: int, avatar_transform: Dictionary = {}) -> int:
    if node is MeshInstance3D:
        var mesh_instance = node as MeshInstance3D
        if mesh_instance.mesh:
            vertex_count = export_mesh_to_obj(mesh_instance.mesh, file, vertex_count, avatar_transform)
    elif node.get_class() == "ImporterMeshInstance3D":
        var importer_mesh = node.get("mesh")
        if importer_mesh:
            var regular_mesh = importer_mesh.get_mesh()
            if regular_mesh:
                vertex_count = export_mesh_to_obj(regular_mesh, file, vertex_count, avatar_transform)
    
    for child in node.get_children():
        vertex_count = export_node_meshes(child, file, vertex_count, avatar_transform)
    
    return vertex_count

func export_mesh_to_obj(mesh: Mesh, file: FileAccess, start_vertex: int, avatar_transform: Dictionary = {}) -> int:
    if not mesh or mesh.get_surface_count() == 0:
        return start_vertex
        
    var arrays = mesh.surface_get_arrays(0)
    if not arrays:
        return start_vertex
        
    var vertices = arrays[Mesh.ARRAY_VERTEX] as PackedVector3Array
    var indices = arrays[Mesh.ARRAY_INDEX] as PackedInt32Array
    
    if not vertices or not indices:
        return start_vertex
    
    for vertex in vertices:
        var blender_coords = Utils.transform_to_blender_coords(vertex)
        
        # Apply same transformation as skeleton
        if avatar_transform.has("scale_factor") and avatar_transform.has("hip_offset"):
            var original_coords = blender_coords
            blender_coords *= avatar_transform.scale_factor
            blender_coords += avatar_transform.hip_offset
            
            if start_vertex == 1:  # Debug first vertex only
                print("DEBUG: Mesh transform - Original: ", original_coords, " -> Transformed: ", blender_coords)
                print("DEBUG: Scale factor: ", avatar_transform.scale_factor, ", Hip offset: ", avatar_transform.hip_offset)
            
        file.store_line("v %f %f %f" % [blender_coords.x, blender_coords.y, blender_coords.z])
    
    for i in range(0, indices.size(), 3):
        var v1 = indices[i] + start_vertex
        var v2 = indices[i+1] + start_vertex
        var v3 = indices[i+2] + start_vertex
        file.store_line("f %d %d %d" % [v1, v3, v2])

    return start_vertex + vertices.size()

func create_mtl_file(mtl_path: String):
    var file = FileAccess.open(mtl_path, FileAccess.WRITE)
    if file:
        file.store_line("# Material file for avatar")
        file.store_line("newmtl default")
        file.store_line("Ka 0.2 0.2 0.2")
        file.store_line("Kd 0.8 0.8 0.8")
        file.store_line("Ks 0.0 0.0 0.0")
        file.close()
