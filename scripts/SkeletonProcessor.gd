class_name SkeletonProcessor
extends RefCounted

const Utils = preload("Utils.gd")
const BoneMapper = preload("BoneMapper.gd")

func extract_armature_skeleton(scene: Node, output_path: String, avatar_transform: Dictionary):
    var skeletons = Utils.find_skeletons(scene)
    if skeletons.size() > 0:
        create_skeleton_obj(skeletons[0], output_path, avatar_transform)
    else:
        print("No skeleton found")

func create_skeleton_obj(skeleton: Skeleton3D, output_path: String, avatar_transform: Dictionary):
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if not file:
        print("Error creating skeleton file")
        return
    
    file.store_line("# Skeleton extracted from GLB")
    var bone_count = skeleton.get_bone_count()
    
    if bone_count == 0:
        file.close()
        print("Empty skeleton found")
        return
    
    print("Found skeleton with ", bone_count, " bones")
    
    if avatar_transform.has("scale_factor") and avatar_transform.scale_factor != 1.0:
        print("Scaling avatar skeleton to match garment proportions")
        var reference_skeleton = load_reference_garment_skeleton()
        if reference_skeleton.vertices.size() > 0:
            write_reference_skeleton_as_target(reference_skeleton, file)
        else:
            print("Could not load reference garment skeleton, falling back to direct conversion")
            write_transformed_avatar_skeleton(skeleton, avatar_transform, file)
    else:
        print("Using direct GLB skeleton conversion")
        write_direct_skeleton_obj(skeleton, file)
    
    file.close()

func load_reference_garment_skeleton():
    var reference_skeleton = {
        "vertices": [],
        "edges": []
    }
    
    var garment_skeleton_path = "garment-data/assets/garments/LCL_Skirt_DressEvening_003/skeleton.obj"
    if FileAccess.file_exists(garment_skeleton_path):
        var file = FileAccess.open(garment_skeleton_path, FileAccess.READ)
        if file:
            while not file.eof_reached():
                var line = file.get_line().strip_edges()
                if line.begins_with("v "):
                    var parts = line.split(" ")
                    if parts.size() >= 4:
                        var vertex = Vector3(float(parts[1]), float(parts[2]), float(parts[3]))
                        reference_skeleton.vertices.append(vertex)
                elif line.begins_with("l "):
                    var parts = line.split(" ")
                    if parts.size() >= 3:
                        var edge = [int(parts[1]) - 1, int(parts[2]) - 1]
                        reference_skeleton.edges.append(edge)
            file.close()
            print("Loaded reference garment skeleton: ", reference_skeleton.vertices.size(), " vertices, ", reference_skeleton.edges.size(), " edges")
    
    return reference_skeleton

func scale_garment_skeleton_to_avatar(reference_skeleton, avatar_skeleton: Skeleton3D):
    print("DEBUG: Scaling garment skeleton to match avatar proportions")
    
    # Find key bones in avatar skeleton for scaling reference
    var bone_mapper = BoneMapper.new()
    var bone_map = bone_mapper.identify_key_bones(avatar_skeleton)
    
    # Calculate avatar scale factors
    var avatar_hip_pos = Utils.calculate_bone_world_position(avatar_skeleton, bone_map.root) if bone_map.root != -1 else Vector3.ZERO
    var avatar_head_pos = Utils.calculate_bone_world_position(avatar_skeleton, bone_map.head) if bone_map.head != -1 else Vector3(0, 2, 0)
    var avatar_height = abs(avatar_head_pos.y - avatar_hip_pos.y)
    
    # Calculate reference garment scale (assume reference is about 1.6m tall)
    var reference_height = 1.6
    var scale_factor = avatar_height / reference_height
    
    print("DEBUG: Avatar height: ", avatar_height, "m, Scale factor: ", scale_factor)
    
    # Scale and position reference skeleton vertices
    var scaled_skeleton = {
        "vertices": [],
        "edges": reference_skeleton.edges  # Keep same connectivity
    }
    
    for vertex in reference_skeleton.vertices:
        # Scale vertex and position relative to avatar hip
        var scaled_vertex = vertex * scale_factor
        # Transform to avatar's coordinate space
        var world_vertex = avatar_hip_pos + scaled_vertex
        var blender_coords = Utils.transform_to_blender_coords(world_vertex)
        scaled_skeleton.vertices.append(blender_coords)
    
    return scaled_skeleton

func write_scaled_garment_skeleton_obj(scaled_skeleton, file):
    print("DEBUG: Writing scaled garment skeleton with ", scaled_skeleton.vertices.size(), " vertices")
    
    # Write vertices
    for vertex in scaled_skeleton.vertices:
        file.store_line("v %f %f %f" % [vertex.x, vertex.y, vertex.z])
    
    # Write edges (connectivity)
    for edge in scaled_skeleton.edges:
        file.store_line("l %d %d" % [edge[0] + 1, edge[1] + 1])

func write_mapped_skeleton_obj(mapped_skeleton, file):
    print("DEBUG: Writing extended skeleton with ", mapped_skeleton.vertices.size(), " vertices")
    
    for vertex in mapped_skeleton.vertices:
        file.store_line("v %f %f %f" % [vertex.x, vertex.y, vertex.z])
    
    for edge in mapped_skeleton.edges:
        file.store_line("l %d %d" % [edge[0] + 1, edge[1] + 1])

func write_direct_skeleton_obj(skeleton: Skeleton3D, file):
    var bone_count = skeleton.get_bone_count()
    
    var world_positions = []
    for i in bone_count:
        var world_pos = Utils.calculate_bone_world_position(skeleton, i)
        var blender_coords = Utils.transform_to_blender_coords(world_pos)
        world_positions.append(blender_coords)
        file.store_line("v %f %f %f" % [blender_coords.x, blender_coords.y, blender_coords.z])
    
    for i in bone_count:
        var parent_idx = skeleton.get_bone_parent(i)
        if parent_idx != -1:
            file.store_line("l %d %d" % [parent_idx + 1, i + 1])

func calculate_avatar_to_garment_transform(scene: Node) -> Dictionary:
    var skeletons = Utils.find_skeletons(scene)
    if skeletons.size() == 0:
        return {"scale_factor": 1.0, "hip_offset": Vector3.ZERO}
    
    var reference_skeleton = load_reference_garment_skeleton()
    if reference_skeleton.vertices.size() == 0:
        return {"scale_factor": 1.0, "hip_offset": Vector3.ZERO}
    
    # Calculate reference garment dimensions
    var ref_min_y = reference_skeleton.vertices[0].y
    var ref_max_y = reference_skeleton.vertices[0].y
    for vertex in reference_skeleton.vertices:
        ref_min_y = min(ref_min_y, vertex.y)
        ref_max_y = max(ref_max_y, vertex.y)
    var garment_height = ref_max_y - ref_min_y
    var garment_hip_y = reference_skeleton.vertices[0].y  # First vertex is hip
    
    # Calculate avatar dimensions  
    var bone_mapper = BoneMapper.new()
    var bone_map = bone_mapper.identify_key_bones(skeletons[0])
    var avatar_hip_pos = Utils.calculate_bone_world_position(skeletons[0], bone_map.root) if bone_map.root != -1 else Vector3.ZERO
    var avatar_head_pos = Utils.calculate_bone_world_position(skeletons[0], bone_map.head) if bone_map.head != -1 else Vector3(0, 2, 0)
    var avatar_height = abs(avatar_head_pos.y - avatar_hip_pos.y)
    
    # Calculate scale factor to match garment skeleton
    var scale_factor = garment_height / avatar_height
    
    # Calculate hip offset to position avatar at garment hip level
    var avatar_hip_blender = Utils.transform_to_blender_coords(avatar_hip_pos)
    var scaled_hip_y = avatar_hip_blender.y * scale_factor
    var hip_offset = Vector3(0, garment_hip_y - scaled_hip_y, 0)
    
    print("DEBUG: Garment height: ", garment_height, "m, Avatar height: ", avatar_height, "m")
    print("DEBUG: Scale factor: ", scale_factor, ", Hip offset: ", hip_offset)
    
    return {"scale_factor": scale_factor, "hip_offset": hip_offset}

func write_reference_skeleton_as_target(reference_skeleton, file):
    print("DEBUG: Writing reference garment skeleton as target (", reference_skeleton.vertices.size(), " vertices)")
    
    # Write the reference skeleton exactly as it is - this is our target
    for vertex in reference_skeleton.vertices:
        file.store_line("v %f %f %f" % [vertex.x, vertex.y, vertex.z])
    
    # Write edges (connectivity)
    for edge in reference_skeleton.edges:
        file.store_line("l %d %d" % [edge[0] + 1, edge[1] + 1])

func write_transformed_avatar_skeleton(skeleton: Skeleton3D, avatar_transform: Dictionary, file):
    print("DEBUG: Writing transformed avatar skeleton")
    var bone_count = skeleton.get_bone_count()
    
    var world_positions = []
    for i in bone_count:
        var world_pos = Utils.calculate_bone_world_position(skeleton, i)
        var blender_coords = Utils.transform_to_blender_coords(world_pos)
        
        # Apply transformation to match garment skeleton
        blender_coords *= avatar_transform.scale_factor
        blender_coords += avatar_transform.hip_offset
        
        world_positions.append(blender_coords)
        file.store_line("v %f %f %f" % [blender_coords.x, blender_coords.y, blender_coords.z])
    
    for i in bone_count:
        var parent_idx = skeleton.get_bone_parent(i)
        if parent_idx != -1:
            file.store_line("l %d %d" % [parent_idx + 1, i + 1])

func create_overlay_skeleton(scene: Node, output_path: String, target_avatar: String, garment_type: String):
    print("Creating overlay skeleton for target avatar: ", target_avatar, " garment type: ", garment_type)
    
    var target_skeleton_path = "garment-data/assets/avatars/" + target_avatar + "/skeleton.obj"
    var target_vertices = []
    var target_edges = []
    
    if load_skeleton_from_file(target_skeleton_path, target_vertices, target_edges):
        print("Loaded target skeleton with ", target_vertices.size(), " vertices and ", target_edges.size(), " edges")
        
        var filtered_edges = filter_edges_for_garment(target_edges, garment_type)
        print("Filtered to ", filtered_edges.size(), " edges for garment type: ", garment_type)
        
        write_overlay_skeleton(output_path, target_vertices, filtered_edges)
    else:
        print("Could not load target skeleton")

func load_skeleton_from_file(file_path: String, vertices: Array, edges: Array) -> bool:
    if not FileAccess.file_exists(file_path):
        print("Target skeleton file not found: ", file_path)
        return false
    
    var file = FileAccess.open(file_path, FileAccess.READ)
    if not file:
        print("Error opening target skeleton file: ", file_path)
        return false
    
    vertices.clear()
    edges.clear()
    
    while not file.eof_reached():
        var line = file.get_line().strip_edges()
        if line.begins_with("v "):
            var parts = line.split(" ")
            if parts.size() >= 4:
                var vertex = Vector3(float(parts[1]), float(parts[2]), float(parts[3]))
                vertices.append(vertex)
        elif line.begins_with("l "):
            var parts = line.split(" ")
            if parts.size() >= 3:
                var edge = [int(parts[1]) - 1, int(parts[2]) - 1]
                edges.append(edge)
    
    file.close()
    return true

func filter_edges_for_garment(edges: Array, garment_type: String) -> Array:
    if garment_type.to_lower() == "skirt":
        return edges
    elif garment_type.to_lower() == "jacket" or garment_type.to_lower() == "shirt":
        return edges
    else:
        return edges

func write_overlay_skeleton(output_path: String, vertices: Array, edges: Array):
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if not file:
        print("Error creating overlay skeleton file: ", output_path)
        return
    
    file.store_line("# Overlay skeleton")
    
    for vertex in vertices:
        file.store_line("v %f %f %f" % [vertex.x, vertex.y, vertex.z])
    
    for edge in edges:
        file.store_line("l %d %d" % [edge[0] + 1, edge[1] + 1])
    
    file.close()
    print("Created overlay skeleton with ", vertices.size(), " vertices and ", edges.size(), " edges")
