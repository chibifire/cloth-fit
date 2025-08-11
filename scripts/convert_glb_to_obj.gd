extends SceneTree

var glb_file: String
var output_name: String
var target_avatar: String = ""
var garment_type: String = ""

func _init():
    var args = OS.get_cmdline_user_args()
    if args.size() < 2:
        print("Usage: godot --headless --script convert_glb_to_obj.gd -- <input.glb> <output_name> [target_avatar] [garment_type]")
        quit(1)
        return
    
    glb_file = args[0]
    output_name = args[1]
    if args.size() > 2:
        target_avatar = args[2]
    if args.size() > 3:
        garment_type = args[3]
    
    call_deferred("_start_conversion")

func _start_conversion():
    print("Loading GLB file: ", glb_file)
    
    var temp_node = Node.new()
    root.add_child(temp_node)
    
    var gltf_document = GLTFDocument.new()
    var gltf_state = GLTFState.new()
    var error = gltf_document.append_from_file(glb_file, gltf_state)
    
    if error != OK:
        print("Error loading GLB file: ", error)
        print("Creating fallback outputs...")
        create_basic_obj(output_name + ".obj")
        create_mtl_file(output_name + ".mtl")
        create_basic_skeleton(output_name + "_skeleton.obj")
        temp_node.queue_free()
        quit(0)
        return
    
    # Generate scene
    var scene = gltf_document.generate_scene(gltf_state)
    
    if not scene:
        print("Error: Could not generate scene from GLB file")
        print("Creating fallback outputs...")
        create_basic_obj(output_name + ".obj")
        create_mtl_file(output_name + ".mtl")
        create_basic_skeleton(output_name + "_skeleton.obj")
        temp_node.queue_free()
        quit(0)
        return
    
    # Add scene to temp node for processing
    temp_node.add_child(scene)
    
    # Export mesh as OBJ
    export_meshes_as_obj(scene, output_name + ".obj")
    create_mtl_file(output_name + ".mtl")
    
    # Extract skeleton and create skeleton.obj
    extract_armature_skeleton(scene, output_name + "_skeleton.obj")
    
    # Create overlay skeleton if target avatar specified
    if target_avatar != "":
        create_overlay_skeleton(scene, output_name + "_" + target_avatar.to_lower() + "_overlay.obj")
    
    # Extract skinning weights and create skin.txt
    extract_skinning_weights(scene, gltf_state, "skin.txt")
    
    print("Converted ", glb_file, " to:")
    print("  - ", output_name, ".obj (mesh)")
    print("  - ", output_name, ".mtl (materials)")
    print("  - ", output_name, "_skeleton.obj (skeleton)")
    if target_avatar != "":
        print("  - ", output_name, "_", target_avatar.to_lower(), "_overlay.obj (overlay skeleton)")
    print("  - skin.txt (skinning weights)")
    
    # Cleanup
    temp_node.queue_free()
    quit(0)

func export_meshes_as_obj(node: Node, output_path: String):
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if not file:
        print("Error creating output file: ", output_path)
        return
    
    file.store_line("# OBJ file exported from Godot")
    file.store_line("mtllib " + output_path.get_basename() + ".mtl")
    var vertex_count = 1
    
    vertex_count = export_node_meshes(node, file, vertex_count)
    file.close()

func export_node_meshes(node: Node, file: FileAccess, vertex_count: int) -> int:
    # Handle both MeshInstance3D and ImporterMeshInstance3D
    if node is MeshInstance3D:
        var mesh_instance = node as MeshInstance3D
        if mesh_instance.mesh:
            vertex_count = export_mesh_to_obj(mesh_instance.mesh, file, vertex_count)
    elif node.get_class() == "ImporterMeshInstance3D":
        # ImporterMeshInstance3D has a different structure
        var importer_mesh = node.get("mesh")
        if importer_mesh:
            # Convert ImporterMesh to regular Mesh
            var regular_mesh = importer_mesh.get_mesh()
            if regular_mesh:
                vertex_count = export_mesh_to_obj(regular_mesh, file, vertex_count)
    
    for child in node.get_children():
        vertex_count = export_node_meshes(child, file, vertex_count)
    
    return vertex_count

func export_mesh_to_obj(mesh: Mesh, file: FileAccess, start_vertex: int) -> int:
    if not mesh or mesh.get_surface_count() == 0:
        return start_vertex
        
    var arrays = mesh.surface_get_arrays(0)
    if not arrays:
        return start_vertex
        
    var vertices = arrays[Mesh.ARRAY_VERTEX] as PackedVector3Array
    var indices = arrays[Mesh.ARRAY_INDEX] as PackedInt32Array
    
    if not vertices or not indices:
        return start_vertex
    
    # Write vertices (match Blender coordinate system: scale and offset to match garment skeleton)
    for vertex in vertices:
        # Transform to match Blender garment skeleton coordinate system
        var blender_x = vertex.x # Map Z to X direction
        var blender_y = -vertex.z         # Map X to Y direction (small scale)
        var blender_z = vertex.y               # Map Y to Z direction
        file.store_line("v %f %f %f" % [blender_x, blender_y, blender_z])
    
    # Write faces (reverse winding order to fix inside-out normals after coordinate transformation)
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

func extract_armature_skeleton(node: Node, output_path: String):
    var skeletons = find_skeletons(node)
    if skeletons.size() > 0:
        create_skeleton_obj(skeletons[0], output_path)
    else:
        print("No skeleton found, creating basic skeleton")
        create_basic_skeleton(output_path)

func find_skeletons(node: Node) -> Array:
    var skeletons = []
    if node is Skeleton3D:
        skeletons.append(node)
    for child in node.get_children():
        skeletons.append_array(find_skeletons(child))
    return skeletons

func create_skeleton_obj(skeleton: Skeleton3D, output_path: String):
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if not file:
        print("Error creating skeleton file, using basic fallback")
        create_basic_skeleton(output_path)
        return
    
    file.store_line("# Skeleton extracted from GLB")
    var bone_count = skeleton.get_bone_count()
    
    if bone_count == 0:
        file.close()
        print("Empty skeleton found, using basic fallback")
        create_basic_skeleton(output_path)
        return
    
    print("Found skeleton with ", bone_count, " bones")
    
    # Load target garment skeleton structure to match
    var target_skeleton = load_target_garment_skeleton()
    if target_skeleton.vertices.size() > 0:
        print("Mapping GLB skeleton to garment skeleton topology...")
        var mapped_skeleton = map_to_garment_skeleton(skeleton, target_skeleton)
        write_mapped_skeleton_obj(mapped_skeleton, file)
    else:
        print("Could not load target skeleton, using direct conversion")
        write_direct_skeleton_obj(skeleton, file)
    
    file.close()

func load_target_garment_skeleton():
    # Load the target garment skeleton structure (LCL_Skirt_DressEvening_003)
    var target_skeleton = {
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
                        target_skeleton.vertices.append(vertex)
                elif line.begins_with("l "):
                    var parts = line.split(" ")
                    if parts.size() >= 3:
                        var edge = [int(parts[1]) - 1, int(parts[2]) - 1]  # Convert to 0-based
                        target_skeleton.edges.append(edge)
            file.close()
            print("Loaded target skeleton: ", target_skeleton.vertices.size(), " vertices, ", target_skeleton.edges.size(), " edges")
    
    return target_skeleton

func map_to_garment_skeleton(glb_skeleton: Skeleton3D, target_skeleton):
    # Create correspondence mapping using the scaling approach
    var mapped_skeleton = {
        "vertices": [],
        "edges": []
    }
    
    # Identify key bones in GLB skeleton
    var bone_map = identify_key_bones(glb_skeleton)
    
    # Calculate scaling factors based on key bone distances
    var scale_factors = calculate_skeleton_scaling(glb_skeleton, target_skeleton, bone_map)
    
    # Map GLB bones to target skeleton structure
    var bone_correspondence = create_bone_correspondence(glb_skeleton, target_skeleton, bone_map, scale_factors)
    
    # Generate the exact 15-vertex structure needed
    for i in range(target_skeleton.vertices.size()):
        if i < bone_correspondence.size() and bone_correspondence[i] != -1:
            var glb_bone_idx = bone_correspondence[i]
            var scaled_pos = scale_bone_position(glb_skeleton, glb_bone_idx, scale_factors)
            mapped_skeleton.vertices.append(scaled_pos)
        else:
            # Use target skeleton vertex if no GLB correspondence
            mapped_skeleton.vertices.append(target_skeleton.vertices[i])
    
    # Copy exact edge structure from target
    mapped_skeleton.edges = target_skeleton.edges.duplicate()
    
    print("Mapped skeleton: ", mapped_skeleton.vertices.size(), " vertices, ", mapped_skeleton.edges.size(), " edges")
    return mapped_skeleton

func identify_key_bones(skeleton: Skeleton3D):
    # Identify key bones based on naming patterns and hierarchy
    var bone_map = {
        "root": -1,
        "spine": -1,
        "chest": -1,
        "neck": -1,
        "head": -1,
        "left_shoulder": -1,
        "right_shoulder": -1,
        "left_arm": -1,
        "right_arm": -1,
        "left_leg": -1,
        "right_leg": -1
    }
    
    for i in range(skeleton.get_bone_count()):
        var bone_name = skeleton.get_bone_name(i).to_lower()
        
        # Root/Hip identification
        if bone_name.contains("hip") or bone_name.contains("pelvis") or bone_name.contains("root"):
            bone_map.root = i
        elif bone_name.contains("spine") or bone_name.contains("torso"):
            bone_map.spine = i
        elif bone_name.contains("chest") or bone_name.contains("upper"):
            bone_map.chest = i
        elif bone_name.contains("neck"):
            bone_map.neck = i
        elif bone_name.contains("head"):
            bone_map.head = i
        elif bone_name.contains("shoulder") and bone_name.contains("l"):
            bone_map.left_shoulder = i
        elif bone_name.contains("shoulder") and bone_name.contains("r"):
            bone_map.right_shoulder = i
        elif bone_name.contains("arm") and bone_name.contains("l"):
            bone_map.left_arm = i
        elif bone_name.contains("arm") and bone_name.contains("r"):
            bone_map.right_arm = i
        elif bone_name.contains("leg") and bone_name.contains("l"):
            bone_map.left_leg = i
        elif bone_name.contains("leg") and bone_name.contains("r"):
            bone_map.right_leg = i
    
    # If root not found, find bone with most connections
    if bone_map.root == -1:
        var max_connections = 0
        for i in range(skeleton.get_bone_count()):
            var connections = 0
            for j in range(skeleton.get_bone_count()):
                if skeleton.get_bone_parent(j) == i:
                    connections += 1
            if connections > max_connections:
                max_connections = connections
                bone_map.root = i
    
    return bone_map

func calculate_skeleton_scaling(glb_skeleton: Skeleton3D, target_skeleton, bone_map):
    # Calculate scaling factors based on key bone distances
    var scale_factors = {
        "global": 1.0,
        "vertical": 1.0,
        "horizontal": 1.0
    }
    
    if bone_map.root != -1 and bone_map.chest != -1:
        # Calculate torso length scaling
        var glb_root_pos = glb_skeleton.get_bone_global_pose(bone_map.root).origin
        var glb_chest_pos = glb_skeleton.get_bone_global_pose(bone_map.chest).origin
        var glb_torso_length = glb_root_pos.distance_to(glb_chest_pos)
        
        # Target skeleton torso length (approximate from vertices 0 to 2)
        if target_skeleton.vertices.size() >= 3:
            var target_torso_length = target_skeleton.vertices[0].distance_to(target_skeleton.vertices[2])
            if glb_torso_length > 0:
                scale_factors.vertical = target_torso_length / glb_torso_length
                scale_factors.global = scale_factors.vertical
    
    # Ensure reasonable scaling bounds
    scale_factors.global = clamp(scale_factors.global, 0.1, 10.0)
    scale_factors.vertical = clamp(scale_factors.vertical, 0.1, 10.0)
    scale_factors.horizontal = scale_factors.global
    
    print("Calculated scale factors: global=", scale_factors.global, " vertical=", scale_factors.vertical)
    return scale_factors

func create_bone_correspondence(glb_skeleton: Skeleton3D, target_skeleton, bone_map, scale_factors):
    # Create correspondence between GLB bones and target skeleton indices
    var correspondence = []
    correspondence.resize(target_skeleton.vertices.size())
    
    # Initialize with -1 (no correspondence)
    for i in range(correspondence.size()):
        correspondence[i] = -1
    
    # Map key bones to specific target indices based on target skeleton structure
    # Target skeleton structure (15 vertices):
    # 0: root, 1: spine_base, 2: spine_top, 3-5: right_arm, 6-8: left_arm, 
    # 9-11: right_leg, 12-14: left_leg
    
    if bone_map.root != -1:
        correspondence[0] = bone_map.root  # Root
        correspondence[1] = bone_map.root  # Spine base (same as root)
    
    if bone_map.chest != -1:
        correspondence[2] = bone_map.chest  # Spine top
    elif bone_map.spine != -1:
        correspondence[2] = bone_map.spine
    
    # Right arm chain (indices 3-5)
    if bone_map.right_shoulder != -1:
        correspondence[3] = bone_map.right_shoulder
    if bone_map.right_arm != -1:
        correspondence[4] = bone_map.right_arm
        correspondence[5] = bone_map.right_arm  # Extend arm
    
    # Left arm chain (indices 6-8)  
    if bone_map.left_shoulder != -1:
        correspondence[6] = bone_map.left_shoulder
    if bone_map.left_arm != -1:
        correspondence[7] = bone_map.left_arm
        correspondence[8] = bone_map.left_arm  # Extend arm
    
    # Right leg chain (indices 9-11)
    if bone_map.right_leg != -1:
        correspondence[9] = bone_map.right_leg
        correspondence[10] = bone_map.right_leg  # Extend leg
        correspondence[11] = bone_map.right_leg  # Foot
    
    # Left leg chain (indices 12-14)
    if bone_map.left_leg != -1:
        correspondence[12] = bone_map.left_leg
        correspondence[13] = bone_map.left_leg  # Extend leg  
        correspondence[14] = bone_map.left_leg  # Foot
    
    return correspondence

func scale_bone_position(glb_skeleton: Skeleton3D, bone_idx: int, scale_factors):
    # Get bone position and scale it to match target coordinate system
    var bone_pose = glb_skeleton.get_bone_global_pose(bone_idx)
    var pos = bone_pose.origin
    
    # Apply scaling
    pos *= scale_factors.global
    
    # Transform to match target coordinate system (around X≈14, Y≈0, Z range)
    var target_x = 14.0 + pos.z * 0.3  # Map Z to X with scaling
    var target_y = pos.x * 0.1          # Map X to Y with small scale
    var target_z = pos.y * 0.5          # Map Y to Z with moderate scale
    
    return Vector3(target_x, target_y, target_z)

func write_mapped_skeleton_obj(mapped_skeleton, file):
    # Write the mapped skeleton with exact target structure
    for vertex in mapped_skeleton.vertices:
        file.store_line("v %f %f %f" % [vertex.x, vertex.y, vertex.z])
    
    # Write edges with exact target connectivity
    for edge in mapped_skeleton.edges:
        file.store_line("l %d %d" % [edge[0] + 1, edge[1] + 1])  # Convert to 1-based

func write_direct_skeleton_obj(skeleton: Skeleton3D, file):
    # Fallback: write skeleton directly with coordinate transformation
    var bone_count = skeleton.get_bone_count()
    
    # Write bone positions as vertices
    for i in bone_count:
        var bone_pose = skeleton.get_bone_global_pose(i)
        var pos = bone_pose.origin
        # Transform to match Blender coordinate system
        var blender_x = 14.0 + pos.z * 0.5
        var blender_y = pos.x * 0.1
        var blender_z = pos.y
        file.store_line("v %f %f %f" % [blender_x, blender_y, blender_z])
    
    # Write bone connections as lines
    for i in bone_count:
        var parent_idx = skeleton.get_bone_parent(i)
        if parent_idx != -1:
            file.store_line("l %d %d" % [parent_idx + 1, i + 1])

func create_overlay_skeleton(scene: Node, output_path: String):
    print("Creating overlay skeleton for target avatar: ", target_avatar, " garment type: ", garment_type)
    
    # Load target skeleton structure
    var target_skeleton_path = "garment-data/assets/avatars/" + target_avatar + "/skeleton.obj"
    var target_vertices = []
    var target_edges = []
    
    if load_skeleton_from_file(target_skeleton_path, target_vertices, target_edges):
        print("Loaded target skeleton with ", target_vertices.size(), " vertices and ", target_edges.size(), " edges")
        
        # Filter edges based on garment type
        var filtered_edges = filter_edges_for_garment(target_edges, garment_type)
        print("Filtered to ", filtered_edges.size(), " edges for garment type: ", garment_type)
        
        # Write overlay skeleton with same vertices but filtered edges
        write_overlay_skeleton(output_path, target_vertices, filtered_edges)
    else:
        print("Could not load target skeleton, creating basic overlay")
        create_basic_skeleton_overlay(output_path)

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
                var edge = [int(parts[1]) - 1, int(parts[2]) - 1]  # Convert from 1-based to 0-based indexing
                edges.append(edge)
    
    file.close()
    return true

func filter_edges_for_garment(edges: Array, garment_type: String) -> Array:
    # For now, return all edges - this can be enhanced with proper filtering logic
    # based on bone names or hierarchy analysis
    
    if garment_type.to_lower() == "skirt":
        # For skirts, we might want to keep only lower body connections
        # This is a simplified implementation - you would enhance this based on actual bone naming
        return edges  # Keep all edges for now
    elif garment_type.to_lower() == "jacket" or garment_type.to_lower() == "shirt":
        # For upper body garments, keep torso/arm connections
        return edges  # Keep all edges for now
    else:
        # Default: keep all edges
        return edges

func write_overlay_skeleton(output_path: String, vertices: Array, edges: Array):
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if not file:
        print("Error creating overlay skeleton file: ", output_path)
        return
    
    file.store_line("# Overlay skeleton for " + target_avatar + " - " + garment_type)
    
    # Write vertices
    for vertex in vertices:
        file.store_line("v %f %f %f" % [vertex.x, vertex.y, vertex.z])
    
    # Write edges
    for edge in edges:
        file.store_line("l %d %d" % [edge[0] + 1, edge[1] + 1])  # Convert back to 1-based indexing
    
    file.close()
    print("Created overlay skeleton with ", vertices.size(), " vertices and ", edges.size(), " edges")

func create_basic_skeleton_overlay(output_path: String):
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if file:
        file.store_line("# Basic overlay skeleton")
        file.store_line("v 0.0 0.0 0.0")
        file.store_line("v 0.0 1.0 0.0")
        file.store_line("v 0.0 1.5 0.0")
        file.store_line("l 1 2")
        file.store_line("l 2 3")
        file.close()

func create_basic_skeleton(output_path: String):
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if file:
        file.store_line("# Basic skeleton")
        file.store_line("v 0.0 0.0 0.0")
        file.store_line("v 0.0 1.0 0.0")
        file.store_line("v 0.0 1.5 0.0")
        file.store_line("l 1 2")
        file.store_line("l 2 3")
        file.close()

func create_basic_obj(output_path: String):
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if file:
        file.store_line("# Basic triangle mesh")
        file.store_line("v -1.0 -1.0 0.0")
        file.store_line("v 1.0 -1.0 0.0")
        file.store_line("v 0.0 1.0 0.0")
        file.store_line("f 1 2 3")
        file.close()

func extract_skinning_weights(scene: Node, gltf_state: GLTFState, output_path: String):
    # Try to extract skinning weights using GLTF state first
    var skins = gltf_state.get_skins()
    var skeletons = find_skeletons(scene)
    
    print("DEBUG: Found ", skins.size(), " skins in GLTF state")
    print("DEBUG: Found ", skeletons.size(), " skeletons")
    
    if skins.size() == 0 or skeletons.size() == 0:
        print("No skins or skeletons foundx")
        return
    
    var skeleton = skeletons[0]
    var bone_count = skeleton.get_bone_count()
    print("DEBUG: Skeleton has ", bone_count, " bones")
    
    # Get the first skin from GLTF state
    var skin = skins[0]
    
    # Try to find mesh instances with actual skinning data
    var mesh_instances = find_mesh_instances(scene)
    var skinned_mesh_data = null
    var vertex_count = 0
    
    # Find the mesh with the most vertices (usually the main character mesh)
    for mesh_instance in mesh_instances:
        var current_mesh = null
        if mesh_instance is MeshInstance3D:
            current_mesh = mesh_instance.mesh
        elif mesh_instance.get_class() == "ImporterMeshInstance3D":
            current_mesh = mesh_instance.get("mesh")
            if current_mesh and current_mesh is ImporterMesh:
                current_mesh = current_mesh.get_mesh()
        
        if current_mesh and current_mesh.get_surface_count() > 0:
            var arrays = current_mesh.surface_get_arrays(0)
            if arrays and arrays[Mesh.ARRAY_VERTEX]:
                var vert_count = arrays[Mesh.ARRAY_VERTEX].size()
                if vert_count > vertex_count:
                    vertex_count = vert_count
                    skinned_mesh_data = arrays
                    print("Found mesh with ", vert_count, " vertices")
    
    if skinned_mesh_data == null or vertex_count == 0:
        print("No suitable mesh found")
        return
    
    # Get skinning data from the mesh arrays
    var weights_array = null
    var bones_array = null
    
    if skinned_mesh_data[Mesh.ARRAY_WEIGHTS] != null:
        weights_array = skinned_mesh_data[Mesh.ARRAY_WEIGHTS]
    if skinned_mesh_data[Mesh.ARRAY_BONES] != null:
        bones_array = skinned_mesh_data[Mesh.ARRAY_BONES]
    
    print("DEBUG: weights_array size: ", weights_array.size() if weights_array else 0)
    print("DEBUG: bones_array size: ", bones_array.size() if bones_array else 0)
    
    if weights_array == null or bones_array == null or weights_array.size() == 0 or bones_array.size() == 0:
        print("No skinning data found in mesh arrays")
        return
    
    # Create bone-centric weight matrix (bone_count x vertex_count)
    var weight_matrix = []
    for bone_idx in range(bone_count):
        var bone_weights = []
        for vertex_idx in range(vertex_count):
            bone_weights.append(0.0)
        weight_matrix.append(bone_weights)
    
    # Fill the weight matrix from vertex-centric data
    # Each vertex has 4 bone influences stored consecutively
    var influences_per_vertex = 4  # Standard for most 3D formats
    
    print("DEBUG: Processing ", vertex_count, " vertices with ", influences_per_vertex, " influences each")
    
    for vertex_idx in range(vertex_count):
        for influence_idx in range(influences_per_vertex):
            var array_idx = vertex_idx * influences_per_vertex + influence_idx
            
            if array_idx < weights_array.size() and array_idx < bones_array.size():
                var weight = weights_array[array_idx]
                var bone_idx = int(bones_array[array_idx])
                
                if weight > 0.0 and bone_idx >= 0 and bone_idx < bone_count:
                    weight_matrix[bone_idx][vertex_idx] = weight
    
    # Write the weight matrix to file
    write_skin_weights_file(weight_matrix, bone_count, vertex_count, output_path)

func find_mesh_instances(node: Node) -> Array:
    var instances = []
    if node is MeshInstance3D or node.get_class() == "ImporterMeshInstance3D":
        instances.append(node)
    for child in node.get_children():
        instances.append_array(find_mesh_instances(child))
    return instances

func write_skin_weights_file(weight_matrix: Array, bone_count: int, vertex_count: int, output_path: String):
    var file = FileAccess.open(output_path, FileAccess.WRITE)
    if not file:
        print("Error creating skin weights file: ", output_path)
        return
    
    file.store_line("# Skin weights extracted from GLB")
    file.store_line("# bone_count: %d" % bone_count)
    file.store_line("# vertex_count: %d" % vertex_count)
    
    # Write weights for each bone
    for bone_idx in range(bone_count):
        var line = ""
        for vertex_idx in range(vertex_count):
            if vertex_idx > 0:
                line += " "
            line += str(weight_matrix[bone_idx][vertex_idx])
        file.store_line(line)
    
    file.close()
    print("Created skin weights file with ", bone_count, " bones and ", vertex_count, " vertices")
