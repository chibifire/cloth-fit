class_name WeightExtractor
extends RefCounted

const Utils = preload("Utils.gd")

func extract_skinning_weights(scene: Node, gltf_state: GLTFState, output_path: String):
    var skins = gltf_state.get_skins()
    var skeletons = Utils.find_skeletons(scene)
    
    print("DEBUG: Found ", skins.size(), " skins in GLTF state")
    print("DEBUG: Found ", skeletons.size(), " skeletons")
    
    if skins.size() == 0 or skeletons.size() == 0:
        print("No skins or skeletons found")
        return
    
    var skeleton = skeletons[0]
    var bone_count = skeleton.get_bone_count()
    print("DEBUG: Skeleton has ", bone_count, " bones")
    
    var mesh_instances = Utils.find_mesh_instances(scene)
    var skinned_mesh_data = find_largest_mesh(mesh_instances)
    
    if not skinned_mesh_data:
        print("No suitable mesh found")
        return
    
    var vertex_count = skinned_mesh_data.vertex_count
    var weights_array = skinned_mesh_data.weights
    var bones_array = skinned_mesh_data.bones
    
    print("DEBUG: weights_array size: ", weights_array.size() if weights_array else 0)
    print("DEBUG: bones_array size: ", bones_array.size() if bones_array else 0)
    
    if not weights_array or not bones_array or weights_array.size() == 0 or bones_array.size() == 0:
        print("No skinning data found in mesh arrays")
        return
    
    var weight_matrix = create_weight_matrix(bone_count, vertex_count, weights_array, bones_array)
    write_skin_weights_file(weight_matrix, bone_count, vertex_count, output_path)

func find_largest_mesh(mesh_instances: Array):
    var largest_mesh_data = null
    var max_vertex_count = 0
    
    for mesh_instance in mesh_instances:
        var current_mesh = extract_mesh_from_instance(mesh_instance)
        
        if current_mesh and current_mesh.get_surface_count() > 0:
            var arrays = current_mesh.surface_get_arrays(0)
            if arrays and arrays[Mesh.ARRAY_VERTEX]:
                var vert_count = arrays[Mesh.ARRAY_VERTEX].size()
                if vert_count > max_vertex_count:
                    max_vertex_count = vert_count
                    largest_mesh_data = {
                        "vertex_count": vert_count,
                        "weights": arrays[Mesh.ARRAY_WEIGHTS],
                        "bones": arrays[Mesh.ARRAY_BONES]
                    }
                    print("Found mesh with ", vert_count, " vertices")
    
    return largest_mesh_data

func extract_mesh_from_instance(mesh_instance):
    if mesh_instance is MeshInstance3D:
        return mesh_instance.mesh
    elif mesh_instance.get_class() == "ImporterMeshInstance3D":
        var importer_mesh = mesh_instance.get("mesh")
        if importer_mesh and importer_mesh is ImporterMesh:
            return importer_mesh.get_mesh()
    return null

func create_weight_matrix(bone_count: int, vertex_count: int, weights_array, bones_array) -> Array:
    var weight_matrix = []
    for bone_idx in range(bone_count):
        var bone_weights = []
        for vertex_idx in range(vertex_count):
            bone_weights.append(0.0)
        weight_matrix.append(bone_weights)
    
    var influences_per_vertex = 4
    
    print("DEBUG: Processing ", vertex_count, " vertices with ", influences_per_vertex, " influences each")
    
    for vertex_idx in range(vertex_count):
        for influence_idx in range(influences_per_vertex):
            var array_idx = vertex_idx * influences_per_vertex + influence_idx
            
            if array_idx < weights_array.size() and array_idx < bones_array.size():
                var weight = weights_array[array_idx]
                var bone_idx = int(bones_array[array_idx])
                
                if weight > 0.0 and bone_idx >= 0 and bone_idx < bone_count:
                    weight_matrix[bone_idx][vertex_idx] = weight
    
    return weight_matrix

func write_skin_weights_file(weight_matrix: Array, bone_count: int, vertex_count: int, output_path: String):
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if not file:
        print("Error creating skin weights file: ", output_path)
        return
    
    file.store_line("# Skin weights extracted from GLB")
    file.store_line("# bone_count: %d" % bone_count)
    file.store_line("# vertex_count: %d" % vertex_count)
    
    for bone_idx in range(bone_count):
        var line = ""
        for vertex_idx in range(vertex_count):
            if vertex_idx > 0:
                line += " "
            line += str(weight_matrix[bone_idx][vertex_idx])
        file.store_line(line)
    
    file.close()
    print("Created skin weights file with ", bone_count, " bones and ", vertex_count, " vertices")
