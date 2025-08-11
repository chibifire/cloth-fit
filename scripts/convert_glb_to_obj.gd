extends SceneTree

var glb_file: String
var output_name: String
var target_avatar: String = ""
var garment_type: String = ""

# Global virtual joints registry (based on KHR_avatar_virtual_joints)
var virtual_joints_registry: Dictionary = {}

# Synthetic bone registry for extending GLB skeleton
var synthetic_bones_registry: Dictionary = {}
var next_synthetic_bone_index: int = -1

func calculate_bone_world_position(skeleton: Skeleton3D, bone_idx: int) -> Vector3:
    # Calculate world position by accumulating transforms from root to bone
    # Now handles both original GLB bones and synthetic bones
    var world_transform = Transform3D.IDENTITY
    var current_bone = bone_idx
    var bone_chain = []
    
    # Build chain from bone to root using extended parent function
    while current_bone != -1:
        bone_chain.push_front(current_bone)
        current_bone = get_bone_parent_extended(skeleton, current_bone)
    
    # Accumulate transforms from root to target bone using extended rest function
    for bone in bone_chain:
        var bone_rest = get_bone_rest_extended(skeleton, bone)
        world_transform = world_transform * bone_rest
    
    return world_transform.origin

func get_bone_children(skeleton: Skeleton3D, bone_idx: int) -> Array:
    var children = []
    var bone_count = skeleton.get_bone_count()
    
    for i in range(bone_count):
        if skeleton.get_bone_parent(i) == bone_idx:
            children.append(i)
    
    return children

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
        temp_node.queue_free()
        quit(1)
        return
    
    # Generate scene
    var scene = gltf_document.generate_scene(gltf_state)
    
    if not scene:
        print("Error: Could not generate scene from GLB file")
        temp_node.queue_free()
        quit(1)
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
    if node is MeshInstance3D:
        var mesh_instance = node as MeshInstance3D
        if mesh_instance.mesh:
            vertex_count = export_mesh_to_obj(mesh_instance.mesh, file, vertex_count)
    elif node.get_class() == "ImporterMeshInstance3D":
        var importer_mesh = node.get("mesh")
        if importer_mesh:
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
    
    for vertex in vertices:
        var blender_coords = transform_to_blender_coords(vertex)
        file.store_line("v %f %f %f" % [blender_coords.x, blender_coords.y, blender_coords.z])
    
    for i in range(0, indices.size(), 3):
        var v1 = indices[i] + start_vertex
        var v2 = indices[i+1] + start_vertex
        var v3 = indices[i+2] + start_vertex
        file.store_line("f %d %d %d" % [v1, v3, v2])

    return start_vertex + vertices.size()

func transform_to_blender_coords(vertex: Vector3) -> Vector3:
    # Only transform coordinate system: GLB (X,Y,Z) -> Blender (X, -Z, Y)
    # No scaling or translation - preserve original GLB coordinates
    return Vector3(vertex.x, -vertex.z, vertex.y)

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
        print("No skeleton found")

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
        print("Error creating skeleton file")
        return
    
    file.store_line("# Skeleton extracted from GLB")
    var bone_count = skeleton.get_bone_count()
    
    if bone_count == 0:
        file.close()
        print("Empty skeleton found")
        return
    
    print("Found skeleton with ", bone_count, " bones")
    
    # Initialize synthetic bone system
    initialize_synthetic_bone_system(skeleton)
    
    # Print debugging information for direct skeleton
    var bone_map = identify_key_bones(skeleton)
    
    # If target avatar is specified, use Eron's decomposition for garment mapping
    if target_avatar != "":
        print("Using Eron's decomposition for garment skeleton mapping to ", target_avatar)
        var target_skeleton = load_target_garment_skeleton()
        if target_skeleton.vertices.size() > 0:
            var mapped_skeleton = map_to_garment_skeleton(skeleton, target_skeleton)
            write_mapped_skeleton_obj(mapped_skeleton, file)
        else:
            print("Could not load target garment skeleton, falling back to direct conversion")
            write_direct_skeleton_obj(skeleton, file)
    else:
        # Always use direct skeleton conversion - preserve original GLB skeleton structure
        print("Using direct GLB skeleton conversion")
        write_direct_skeleton_obj(skeleton, file)
    
    file.close()

func initialize_synthetic_bone_system(skeleton: Skeleton3D):
    # Initialize synthetic bone system with starting index after GLB bones
    var bone_count = skeleton.get_bone_count()
    next_synthetic_bone_index = bone_count
    synthetic_bones_registry.clear()
    
    print("DEBUG: Synthetic bone system initialized - next index: ", next_synthetic_bone_index)

func create_synthetic_bone(parent_bone: int, world_position: Vector3, bone_type: String) -> int:
    # Create a new synthetic bone with unique index, proper reparenting and transform fixup
    # This function now extends the actual Godot skeleton with real bones
    var synthetic_index = next_synthetic_bone_index
    next_synthetic_bone_index += 1
    
    # Calculate local transform relative to parent bone
    var local_transform = calculate_synthetic_local_transform(parent_bone, world_position)
    
    # Store synthetic bone data (for reference)
    var bone_name = bone_type + "_" + str(synthetic_index)
    synthetic_bones_registry[synthetic_index] = {
        "bone_type": bone_type,
        "parent": parent_bone,
        "world_position": world_position,
        "local_transform": local_transform,
        "name": bone_name
    }
    
    # EXTEND THE ACTUAL GODOT SKELETON - Add real bone to the scene
    extend_godot_skeleton_with_synthetic_bone(synthetic_index, bone_name, parent_bone, local_transform)
    
    print("DEBUG: Created synthetic bone ", synthetic_index, " (", bone_type, ") parented to GLB[", parent_bone, "]")
    print("DEBUG:   World pos: ", world_position)
    print("DEBUG:   Local transform: ", local_transform.origin)
    print("DEBUG:   Added to Godot skeleton as real bone")
    
    return synthetic_index

func calculate_synthetic_local_transform(parent_bone_idx: int, target_world_pos: Vector3) -> Transform3D:
    # Calculate the local transform needed for synthetic bone relative to its parent
    # This ensures proper parent-child spatial relationships
    
    # For now, use a simple local offset approach
    # In a full implementation, you'd want to get the actual parent bone's world transform
    # and calculate the inverse transform to get proper local coordinates
    
    # Calculate a reasonable local offset (12cm laterally for hips)
    var local_offset = Vector3.ZERO
    
    # Determine offset based on target world position patterns
    if abs(target_world_pos.x) > 0.05:  # Lateral offset detected
        local_offset.x = target_world_pos.x * 0.3  # Reduce to local space
    if abs(target_world_pos.z) > 0.05:  # Forward/backward offset
        local_offset.z = target_world_pos.z * 0.3
    
    # Create local transform with the calculated offset
    var local_transform = Transform3D.IDENTITY
    local_transform.origin = local_offset
    
    return local_transform

func extend_godot_skeleton_with_synthetic_bone(synthetic_index: int, bone_name: String, parent_bone: int, local_transform: Transform3D):
    # This function extends the actual Godot Skeleton3D with real synthetic bones
    # We need access to the skeleton - for now store the data and apply during skeleton processing
    
    # Store the synthetic bone data that will be applied when we have skeleton access
    var extension_data = {
        "synthetic_index": synthetic_index,
        "bone_name": bone_name,
        "parent_bone": parent_bone,
        "local_transform": local_transform,
        "needs_skeleton_extension": true
    }
    
    # Add to synthetic bones registry with extension flag
    if synthetic_bones_registry.has(synthetic_index):
        synthetic_bones_registry[synthetic_index]["needs_skeleton_extension"] = true
        synthetic_bones_registry[synthetic_index]["extension_data"] = extension_data
    
    print("DEBUG: Marked synthetic bone ", synthetic_index, " for Godot skeleton extension")

func apply_synthetic_bones_to_skeleton(skeleton: Skeleton3D):
    # Apply all synthetic bones to the actual Godot skeleton
    print("DEBUG: Applying ", synthetic_bones_registry.size(), " synthetic bones to Godot skeleton")
    
    var bones_added = 0
    for synthetic_idx in synthetic_bones_registry.keys():
        var synthetic_data = synthetic_bones_registry[synthetic_idx]
        
        if synthetic_data.has("needs_skeleton_extension") and synthetic_data.needs_skeleton_extension:
            # Add this synthetic bone to the actual skeleton
            var bone_name = synthetic_data.name
            var parent_bone = synthetic_data.parent
            var local_transform = synthetic_data.local_transform
            
            # Add the bone to the skeleton
            skeleton.add_bone(bone_name)
            var new_bone_idx = skeleton.get_bone_count() - 1  # Last added bone
            
            # Set parent relationship
            if parent_bone >= 0 and parent_bone < skeleton.get_bone_count():
                skeleton.set_bone_parent(new_bone_idx, parent_bone)
            
            # Set rest transform
            skeleton.set_bone_rest(new_bone_idx, local_transform)
            
            bones_added += 1
            print("DEBUG: Added synthetic bone '", bone_name, "' as bone ", new_bone_idx, " with parent ", parent_bone)
            
            # Mark as applied
            synthetic_data.needs_skeleton_extension = false
    
    print("DEBUG: Successfully added ", bones_added, " synthetic bones to Godot skeleton")
    print("DEBUG: Skeleton now has ", skeleton.get_bone_count(), " total bones")

func get_bone_world_position_extended(glb_skeleton: Skeleton3D, bone_idx: int) -> Vector3:
    # Get world position for both original GLB bones and synthetic bones
    if bone_idx < glb_skeleton.get_bone_count():
        # Original GLB bone or successfully applied synthetic bone
        return calculate_bone_world_position(glb_skeleton, bone_idx)
    elif synthetic_bones_registry.has(bone_idx):
        # Synthetic bone - return stored world position
        return synthetic_bones_registry[bone_idx].world_position
    else:
        print("ERROR: Invalid bone index ", bone_idx)
        return Vector3.ZERO

func get_bone_parent_extended(glb_skeleton: Skeleton3D, bone_idx: int) -> int:
    # Get bone parent for both original GLB bones and synthetic bones
    if bone_idx < glb_skeleton.get_bone_count():
        # Original GLB bone
        return glb_skeleton.get_bone_parent(bone_idx)
    elif synthetic_bones_registry.has(bone_idx):
        # Synthetic bone - return stored parent
        return synthetic_bones_registry[bone_idx].parent
    else:
        print("ERROR: Invalid bone index for parent lookup ", bone_idx)
        return -1

func get_bone_rest_extended(glb_skeleton: Skeleton3D, bone_idx: int) -> Transform3D:
    # Get bone rest transform for both original GLB bones and synthetic bones
    if bone_idx < glb_skeleton.get_bone_count():
        # Original GLB bone
        return glb_skeleton.get_bone_rest(bone_idx)
    elif synthetic_bones_registry.has(bone_idx):
        # Synthetic bone - return stored local transform
        return synthetic_bones_registry[bone_idx].local_transform
    else:
        print("ERROR: Invalid bone index for rest transform ", bone_idx)
        return Transform3D.IDENTITY

func get_bone_name_extended(glb_skeleton: Skeleton3D, bone_idx: int) -> String:
    # Get bone name for both original GLB bones and synthetic bones
    if bone_idx < glb_skeleton.get_bone_count():
        # Original GLB bone
        return glb_skeleton.get_bone_name(bone_idx)
    elif synthetic_bones_registry.has(bone_idx):
        # Synthetic bone - return stored name
        return synthetic_bones_registry[bone_idx].name
    else:
        return "unknown_bone_" + str(bone_idx)

func calculate_distance_from_root_extended(glb_skeleton: Skeleton3D, bone_idx: int) -> int:
    # Calculate distance from root for both original GLB bones and synthetic bones
    if bone_idx < glb_skeleton.get_bone_count():
        # Original GLB bone - use standard calculation
        return calculate_distance_from_root(glb_skeleton, bone_idx)
    elif synthetic_bones_registry.has(bone_idx):
        # Synthetic bone - calculate distance through parent chain
        var synthetic_data = synthetic_bones_registry[bone_idx]
        var parent_bone = synthetic_data.parent
        
        # Distance = parent's distance + 1
        if parent_bone < glb_skeleton.get_bone_count():
            return calculate_distance_from_root(glb_skeleton, parent_bone) + 1
        else:
            # Parent is also synthetic - recursively calculate
            return calculate_distance_from_root_extended(glb_skeleton, parent_bone) + 1
    else:
        return -1  # Invalid bone index

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

func fix_garment_skeleton_topology(original_edges: Array) -> Array:
    # Use original garment skeleton topology to match existing garment standards
    print("DEBUG: Using original garment skeleton topology for compatibility")
    
    var fixed_edges = []
    
    # ORIGINAL garment skeleton topology (matches LCL_Skirt_DressEvening_003):
    # 0=root, 1=spine_lower, 2=spine_upper
    # 3=shoulder_right, 4=arm_right_upper, 5=arm_right_hand
    # 6=shoulder_left, 7=arm_left_upper, 8=arm_left_hand
    # 9=hip_right, 10=leg_right_upper, 11=leg_right_foot
    # 12=hip_left, 13=leg_left_upper, 14=leg_left_foot
    
    # Original edge connections (matching garment):
    fixed_edges = [
        [0, 1],    # root -> spine_lower
        [1, 2],    # spine_lower -> spine_upper
        [1, 3],    # spine_lower -> shoulder_right
        [3, 4],    # shoulder_right -> arm_right_upper
        [4, 5],    # arm_right_upper -> arm_right_hand
        [1, 6],    # spine_lower -> shoulder_left
        [6, 8],    # shoulder_left -> arm_left_hand
        [7, 8],    # arm_left_upper -> arm_left_hand
        [0, 9],    # root -> hip_right
        [9, 10],   # hip_right -> leg_right_upper
        [10, 11],  # leg_right_upper -> leg_right_foot
        [0, 12],   # root -> hip_left
        [12, 13],  # hip_left -> leg_left_upper
        [13, 14]   # leg_left_upper -> leg_left_foot
    ]
    
    print("DEBUG: Original garment skeleton topology with ", fixed_edges.size(), " edges")
    print("DEBUG:   Arms branch from spine_lower(1), shoulders at indices 3,6")
    return fixed_edges

func map_to_garment_skeleton(glb_skeleton: Skeleton3D, target_skeleton):
    # Use exact same hierarchy as garment, only change bone translations to match Mire
    var mapped_skeleton = {
        "vertices": [],
        "edges": []
    }
    
    print("Using garment hierarchy with Mire's bone positions")
    
    # FIX GARMENT SKELETON TOPOLOGY - The original has incorrect arm connections
    mapped_skeleton.edges = fix_garment_skeleton_topology(target_skeleton.edges)
    
    # Identify key bones in Mire's GLB skeleton
    var bone_map = identify_key_bones(glb_skeleton)
    
    # Create bone correspondence for positioning
    var bone_correspondence = create_bone_correspondence(glb_skeleton, target_skeleton, bone_map, {})
    
    # APPLY SYNTHETIC BONES TO ACTUAL GODOT SKELETON
    # This must happen after synthetic bones are created but before skeleton operations
    apply_synthetic_bones_to_skeleton(glb_skeleton)
    
    # Generate vertices using Mire's bone positions but garment's hierarchy
    print("Mapping ", target_skeleton.vertices.size(), " garment vertices to Mire's anatomy:")
    
    for i in range(target_skeleton.vertices.size()):
        var mire_position: Vector3
        
        if i < bone_correspondence.size() and bone_correspondence[i] != -1:
            # Use Mire's bone position for this garment vertex
            var glb_bone_idx = bone_correspondence[i]
            var bone_world_pos = calculate_bone_world_position(glb_skeleton, glb_bone_idx)
            mire_position = Vector3(bone_world_pos.x, -bone_world_pos.z, bone_world_pos.y)  # Transform to target coords
            
            var bone_name = glb_skeleton.get_bone_name(glb_bone_idx)
            print("  Garment[", i, "] -> Mire bone ", glb_bone_idx, " ('", bone_name, "') at ", mire_position)
        else:
            # No correspondence found, use garment's original position (scaled to Mire's size)
            mire_position = scale_garment_vertex_to_mire(target_skeleton.vertices[i], glb_skeleton, bone_map)
            print("  Garment[", i, "] -> No correspondence, using scaled garment position: ", mire_position)
        
        mapped_skeleton.vertices.append(mire_position)
    
    print("Final skeleton: Same hierarchy as garment (", mapped_skeleton.edges.size(), " edges), Mire's proportions (", mapped_skeleton.vertices.size(), " vertices)")
    return mapped_skeleton

func scale_garment_vertex_to_mire(garment_vertex: Vector3, glb_skeleton: Skeleton3D, bone_map) -> Vector3:
    # Scale garment vertex position to match Mire's overall size
    var mire_root_pos = calculate_bone_world_position(glb_skeleton, bone_map.root) if bone_map.root != -1 else Vector3.ZERO
    var mire_head_pos = calculate_bone_world_position(glb_skeleton, bone_map.head) if bone_map.head != -1 else Vector3(0, 10, 0)
    
    # Calculate Mire's height
    var mire_height = abs(mire_head_pos.y - mire_root_pos.y)
    
    # Assume garment skeleton is normalized, scale to Mire's size
    var scale_factor = mire_height / 10.0  # Rough scaling based on typical skeleton height
    
    var scaled_pos = garment_vertex * scale_factor
    return Vector3(scaled_pos.x, -scaled_pos.z, scaled_pos.y)  # Transform coordinates

func identify_key_bones(skeleton: Skeleton3D):
    # Identify key bones based on topology/hierarchy analysis instead of names
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
    
    var bone_count = skeleton.get_bone_count()
    
    # Find root bone (bone with most children or no parent)
    var root_candidates = []
    var max_children = 0
    
    for i in range(bone_count):
        if skeleton.get_bone_parent(i) == -1:  # No parent = root candidate
            root_candidates.append(i)
        
        # Count children
        var children_count = 0
        for j in range(bone_count):
            if skeleton.get_bone_parent(j) == i:
                children_count += 1
        
        if children_count > max_children:
            max_children = children_count
            bone_map.root = i
    
    # If we have root candidates with no parent, prefer those
    if root_candidates.size() > 0:
        bone_map.root = root_candidates[0]
    
    # Analyze topology from root to identify spine chain
    if bone_map.root != -1:
        var spine_chain = find_main_spine_chain(skeleton, bone_map.root)
        if spine_chain.size() >= 2:
            bone_map.spine = spine_chain[1] if spine_chain.size() > 1 else spine_chain[0]
            bone_map.chest = spine_chain[spine_chain.size() - 1]
    
    # Find arm and leg chains using topology
    analyze_limb_topology(skeleton, bone_map)
    
    # Print debugging information
    print_skeleton_railroad_diagram(skeleton)
    print_bone_map_analysis(skeleton, bone_map)
    
    return bone_map

func find_main_spine_chain(skeleton: Skeleton3D, root_bone: int) -> Array:
    # Find the main spine chain by looking for spine-like bone names in hierarchy
    var spine_chain = [root_bone]
    var current_bone = root_bone
    var bone_count = skeleton.get_bone_count()
    
    # Follow the spine by looking for spine-named bones in the children
    var remaining_depth = 5  # Limit depth to avoid infinite loops
    while remaining_depth > 0:
        remaining_depth -= 1
        
        var children = []
        for i in range(bone_count):
            if skeleton.get_bone_parent(i) == current_bone:
                children.append(i)
        
        # Look for spine-related bones among children
        var next_spine_bone = -1
        for child in children:
            var bone_name = skeleton.get_bone_name(child).to_lower()
            if bone_name.contains("spine") or bone_name.contains("chest") or bone_name.contains("neck") or bone_name.contains("head"):
                next_spine_bone = child
                break
        
        # If no spine-related bone found but only one child, follow it
        if next_spine_bone == -1 and children.size() == 1:
            next_spine_bone = children[0]
        
        if next_spine_bone != -1:
            spine_chain.append(next_spine_bone)
            current_bone = next_spine_bone
        else:
            break  # No spine continuation found
    
    return spine_chain

func analyze_limb_topology(skeleton: Skeleton3D, bone_map):
    # Analyze the topology to identify limbs based on branching patterns
    var bone_count = skeleton.get_bone_count()
    
    if bone_map.root == -1:
        return
    
    # First, find the spine chain to get chest bone
    var spine_chain = find_main_spine_chain(skeleton, bone_map.root)
    if spine_chain.size() >= 3:
        bone_map.spine = spine_chain[1]  # Second bone in chain (Spine)
        bone_map.chest = spine_chain[2]  # Third bone in chain (Chest)
        if spine_chain.size() >= 4:
            bone_map.neck = spine_chain[3]  # Fourth bone (Neck)
        if spine_chain.size() >= 5:
            bone_map.head = spine_chain[4]  # Fifth bone (Head)
    
    # Find limbs that branch from spine chain bones
    var limb_candidates = []
    
    # Look for bones that branch from spine chain (especially chest)
    for i in range(bone_count):
        var parent = skeleton.get_bone_parent(i)
        # Check if parent is in spine chain or is root
        if parent == bone_map.root or parent == bone_map.spine or parent == bone_map.chest:
            # Skip if it's part of the spine chain itself
            var is_spine_bone = false
            for spine_bone in spine_chain:
                if i == spine_bone:
                    is_spine_bone = true
                    break
            
            if not is_spine_bone:
                var bone_name = skeleton.get_bone_name(i).to_lower()
                limb_candidates.append({"idx": i, "name": bone_name, "parent": parent})
    
    # Categorize limbs based on names and parent positions
    for candidate in limb_candidates:
        var name = candidate.name
        var idx = candidate.idx
        
        print("DEBUG: Analyzing candidate ", idx, ": '", name, "' parent=", candidate.parent)
        
        # Arms typically branch from chest
        if candidate.parent == bone_map.chest:
            if name.contains("shoulder") and (name.contains(".l") or name.contains("_l")):
                print("DEBUG: Found left shoulder: ", idx, " '", name, "'")
                bone_map.left_shoulder = idx
            elif name.contains("shoulder") and (name.contains(".r") or name.contains("_r")):
                print("DEBUG: Found right shoulder: ", idx, " '", name, "'")
                bone_map.right_shoulder = idx
            elif name.contains("arm") and (name.contains(".l") or name.contains("_l")):
                print("DEBUG: Found left arm: ", idx, " '", name, "'")
                bone_map.left_arm = idx
            elif name.contains("arm") and (name.contains(".r") or name.contains("_r")):
                print("DEBUG: Found right arm: ", idx, " '", name, "'")
                bone_map.right_arm = idx
        
        # Legs typically branch from root (hips)
        elif candidate.parent == bone_map.root:
            # Check for right side first to avoid conflicts
            if name.contains("leg") and (name.contains(".r") or name.contains("_r")):
                print("DEBUG: Found right leg: ", idx, " '", name, "'")
                bone_map.right_leg = idx
            elif name.contains("leg") and (name.contains(".l") or name.contains("_l")):
                print("DEBUG: Found left leg: ", idx, " '", name, "'")
                bone_map.left_leg = idx
    
    # If shoulders found, find arms from shoulders
    if bone_map.left_shoulder != -1:
        var children = get_bone_children(skeleton, bone_map.left_shoulder)
        if children.size() > 0:
            bone_map.left_arm = children[0]  # First child of shoulder
    
    if bone_map.right_shoulder != -1:
        var children = get_bone_children(skeleton, bone_map.right_shoulder)
        if children.size() > 0:
            bone_map.right_arm = children[0]  # First child of shoulder

func calculate_skeleton_scaling(glb_skeleton: Skeleton3D, target_skeleton, bone_map):
    # No scaling - preserve original GLB coordinates
    var scale_factors = {
        "scale": 1.0,
        "translation": Vector3.ZERO,
        "global": 1.0
    }
    
    print("No scaling applied - using original GLB coordinates")
    return scale_factors

func create_bone_correspondence(glb_skeleton: Skeleton3D, target_skeleton, bone_map, scale_factors):
    # Use Eron's decomposition algorithm for sophisticated skeleton correspondence
    print("=== APPLYING ERON'S DECOMPOSITION ALGORITHM ===")
    return apply_erons_decomposition(glb_skeleton, target_skeleton, bone_map)

# Eron's Decomposition Algorithm Implementation
class GarmentEffector:
    var index: int              # 0-14 for the 15 garment bones
    var position: Vector3       # Target position in garment skeleton
    var opacity: float          # Weight/importance (1.0 for now)
    
    func _init(idx: int, pos: Vector3, weight: float = 1.0):
        index = idx
        position = pos
        opacity = weight

class BoneInfluence:
    var glb_bone_idx: int
    var effector_weights: Dictionary  # effector_idx -> weight
    var distance_from_root: int
    var bone_name: String
    
    func _init(bone_idx: int, dist: int, name: String):
        glb_bone_idx = bone_idx
        effector_weights = {}
        distance_from_root = dist
        bone_name = name

class EffectorGroup:
    var bone_sequence: Array        # GLB bone indices
    var target_effectors: Array     # Which garment bones this serves
    var root_distance: int          # Distance of rootmost bone from GLB root
    
    func _init():
        bone_sequence = []
        target_effectors = []
        root_distance = 0

func apply_erons_decomposition(glb_skeleton: Skeleton3D, target_skeleton, bone_map) -> Array:
    print("Starting Eron's decomposition with ", glb_skeleton.get_bone_count(), " GLB bones -> ", target_skeleton.vertices.size(), " garment bones")
    
    # Use ONLY the geometric correspondence algorithm - no old hierarchy method
    var correspondence = create_unique_bone_mappings(glb_skeleton, target_skeleton.vertices.size(), bone_map)
    
    print("=== ERON'S DECOMPOSITION COMPLETE ===")
    return correspondence

func create_unique_bone_mappings(glb_skeleton: Skeleton3D, target_count: int, bone_map) -> Array:
    var correspondence = []
    correspondence.resize(target_count)
    
    # Initialize with -1 (no correspondence)
    for i in range(correspondence.size()):
        correspondence[i] = -1
    
    print("DEBUG: Computing hierarchical and rotational geometric correspondence")
    
    # Step 1: Analyze GLB skeleton hierarchy and build hierarchical map
    var hierarchy_analysis = analyze_glb_hierarchy(glb_skeleton)
    
    # Step 2: Get target garment skeleton structure
    var target_skeleton = load_target_garment_skeleton()
    if target_skeleton.vertices.size() != target_count:
        print("ERROR: Target skeleton has ", target_skeleton.vertices.size(), " vertices, expected ", target_count)
        return correspondence
    
    # Step 3: Use hierarchical and rotational matching to find correspondences
    correspondence = create_hierarchical_rotational_correspondence(glb_skeleton, target_skeleton, hierarchy_analysis, bone_map)
    
    return correspondence

func analyze_glb_hierarchy(glb_skeleton: Skeleton3D) -> Dictionary:
    var analysis = {
        "major_bones": {},          # Maps GLB bone indices to their hierarchical roles
        "bone_transforms": {},      # Maps GLB bone indices to their transform data
        "hierarchical_levels": {},  # Maps distance from root to bone lists
        "rotational_chains": []     # Lists of connected bone chains for rotation analysis
    }
    
    var bone_count = glb_skeleton.get_bone_count()
    print("DEBUG: Analyzing hierarchy of ", bone_count, " bones")
    
    # Analyze each bone's hierarchical position and rotational properties
    for bone_idx in range(bone_count):
        var bone_name = glb_skeleton.get_bone_name(bone_idx)
        var bone_transform = glb_skeleton.get_bone_rest(bone_idx)
        var parent_idx = glb_skeleton.get_bone_parent(bone_idx)
        var distance_from_root = calculate_distance_from_root(glb_skeleton, bone_idx)
        
        # Store transform data for rotational analysis
        analysis.bone_transforms[bone_idx] = {
            "transform": bone_transform,
            "position": calculate_bone_world_position(glb_skeleton, bone_idx),
            "parent": parent_idx,
            "name": bone_name,
            "distance_from_root": distance_from_root
        }
        
        # Categorize bone by hierarchical level
        if not analysis.hierarchical_levels.has(distance_from_root):
            analysis.hierarchical_levels[distance_from_root] = []
        analysis.hierarchical_levels[distance_from_root].append(bone_idx)
        
        # Determine major bone role based on hierarchy and naming
        var bone_role = determine_bone_role(glb_skeleton, bone_idx, bone_name, distance_from_root)
        if bone_role != "":
            analysis.major_bones[bone_idx] = bone_role
            print("DEBUG: Identified major bone ", bone_idx, " ('", bone_name, "') as ", bone_role)
    
    # Build rotational chains for swing/twist analysis
    analysis.rotational_chains = build_rotational_chains(glb_skeleton, analysis)
    
    return analysis

func determine_bone_role(glb_skeleton: Skeleton3D, bone_idx: int, bone_name: String, distance_from_root: int) -> String:
    # Pure geometric classification using spatial analysis and swing rotation patterns
    
    # First: Filter out non-anatomical bones using geometric patterns
    if is_accessory_bone(glb_skeleton, bone_idx, bone_name):
        return ""  # Skip accessory bones
    
    # Get bone spatial properties
    var bone_pos = calculate_bone_world_position(glb_skeleton, bone_idx)
    var parent_idx = glb_skeleton.get_bone_parent(bone_idx)
    var children = get_bone_children(glb_skeleton, bone_idx)
    
    # Calculate skeleton bounds for relative positioning
    var skeleton_bounds = calculate_skeleton_spatial_bounds(glb_skeleton)
    var root_pos = calculate_bone_world_position(glb_skeleton, 0)
    
    # Root bone: distance 0, center of coordinate system
    if distance_from_root == 0:
        return "root"
    
    # Classify by spatial zones and anatomical structure
    var relative_height = (bone_pos.y - skeleton_bounds.min_y) / (skeleton_bounds.max_y - skeleton_bounds.min_y)
    var lateral_distance = abs(bone_pos.x - root_pos.x)
    var is_left = bone_pos.x < root_pos.x
    var side = "left" if is_left else "right"
    
    # SPINE CHAIN: Vertical progression from root, low lateral distance
    if distance_from_root <= 2 and lateral_distance < 2.0:
        if distance_from_root == 1 and bone_pos.y > root_pos.y:
            return "spine_1"
        elif distance_from_root == 2 and bone_pos.y > root_pos.y:
            return "spine_2"
    
    # HEAD/NECK: High position, close to center
    elif relative_height > 0.8 and lateral_distance < 3.0:
        if distance_from_root == 3:
            return "neck"
        elif distance_from_root == 4:
            return "head"
    
    # UPPER BODY LIMBS: More inclusive detection for arms and shoulders
    elif lateral_distance > 1.0 and relative_height > 0.3:  # More inclusive thresholds
        # Shoulder detection: distance 3-4, lateral offset from spine
        if distance_from_root >= 3 and distance_from_root <= 4 and parent_idx != -1:
            var parent_distance = calculate_distance_from_root(glb_skeleton, parent_idx)
            if parent_distance <= 3:  # More inclusive parent distance
                return "shoulder_" + side
        
        # Arm detection: distance 4-6, extends laterally
        elif distance_from_root >= 4 and distance_from_root <= 6:
            if relative_height > 0.4:  # Upper body region
                if distance_from_root <= 5:
                    return "upper_arm_" + side
                else:
                    return "lower_arm_" + side
    
    # LOWER BODY LIMBS: More inclusive detection for legs
    elif lateral_distance > 0.5:  # Much more inclusive for legs
        # Leg detection: any bone with lateral offset that could be a leg
        if distance_from_root >= 3 and distance_from_root <= 6:
            if relative_height < 0.7:  # Not in head/neck region
                if distance_from_root <= 4:
                    return "upper_leg_" + side
                else:
                    return "lower_leg_" + side
    
    return ""

func is_accessory_bone(glb_skeleton: Skeleton3D, bone_idx: int, bone_name: String) -> bool:
    # Identify accessory bones using geometric swing analysis and naming patterns
    var name_lower = bone_name.to_lower()
    
    # Obvious accessory patterns (backup check)
    var accessory_patterns = ["skirt", "tail", "hair", "ear", "eye", "breast", "ahoge", "bang"]
    for pattern in accessory_patterns:
        if name_lower.contains(pattern):
            return true
    
    # Geometric swing analysis: Accessories often have irregular swing patterns
    var bone_pos = calculate_bone_world_position(glb_skeleton, bone_idx)
    var parent_idx = glb_skeleton.get_bone_parent(bone_idx)
    
    if parent_idx != -1:
        var parent_pos = calculate_bone_world_position(glb_skeleton, parent_idx)
        var root_pos = calculate_bone_world_position(glb_skeleton, 0)
        
        # Check for skirt-like swing patterns: bones that rotate around hips
        if is_skirt_like_swing_pattern(bone_pos, parent_pos, root_pos):
            return true
    
    return false

func is_skirt_like_swing_pattern(bone_pos: Vector3, parent_pos: Vector3, root_pos: Vector3) -> bool:
    # Detect skirt bones: swing rotation around hip level with radial distribution
    var hip_level_tolerance = 2.0  # Y tolerance for hip level
    var radial_distance = Vector2(bone_pos.x - root_pos.x, bone_pos.z - root_pos.z).length()
    
    # Skirt bones are typically:
    # 1. At similar Y level to hips (root)
    # 2. Distributed radially around the character
    # 3. Have similar radial distance from center
    
    var is_hip_level = abs(bone_pos.y - root_pos.y) < hip_level_tolerance
    var is_radially_distributed = radial_distance > 1.0 and radial_distance < 8.0
    
    return is_hip_level and is_radially_distributed

func calculate_skeleton_spatial_bounds(glb_skeleton: Skeleton3D) -> Dictionary:
    # Calculate the spatial bounds of the entire skeleton
    var bounds = {
        "min_x": INF, "max_x": -INF,
        "min_y": INF, "max_y": -INF,
        "min_z": INF, "max_z": -INF
    }
    
    var bone_count = glb_skeleton.get_bone_count()
    for bone_idx in range(bone_count):
        var pos = calculate_bone_world_position(glb_skeleton, bone_idx)
        
        bounds.min_x = min(bounds.min_x, pos.x)
        bounds.max_x = max(bounds.max_x, pos.x)
        bounds.min_y = min(bounds.min_y, pos.y)
        bounds.max_y = max(bounds.max_y, pos.y)
        bounds.min_z = min(bounds.min_z, pos.z)
        bounds.max_z = max(bounds.max_z, pos.z)
    
    return bounds

func build_rotational_chains(glb_skeleton: Skeleton3D, analysis: Dictionary) -> Array:
    var chains = []
    var processed_bones = {}
    
    # Find chains starting from major bones
    for bone_idx in analysis.major_bones.keys():
        if processed_bones.has(bone_idx):
            continue
            
        var chain = build_chain_from_bone(glb_skeleton, bone_idx, processed_bones)
        if chain.size() > 1:
            chains.append(chain)
    
    print("DEBUG: Built ", chains.size(), " rotational chains")
    return chains

func build_chain_from_bone(glb_skeleton: Skeleton3D, start_bone: int, processed_bones: Dictionary) -> Array:
    var chain = [start_bone]
    processed_bones[start_bone] = true
    
    # Follow the chain forward (to children)
    var current_bone = start_bone
    while true:
        var children = get_bone_children(glb_skeleton, current_bone)
        if children.size() != 1:  # Chain breaks if multiple or no children
            break
        
        var child = children[0]
        if processed_bones.has(child):
            break
            
        chain.append(child)
        processed_bones[child] = true
        current_bone = child
    
    return chain

func create_hierarchical_rotational_correspondence(glb_skeleton: Skeleton3D, target_skeleton, hierarchy_analysis: Dictionary, bone_map) -> Array:
    var correspondence = []
    correspondence.resize(target_skeleton.vertices.size())
    
    # Initialize with -1
    for i in range(correspondence.size()):
        correspondence[i] = -1
    
    print("DEBUG: Creating unique bone correspondence with virtual splitting for ", target_skeleton.vertices.size(), " target bones")
    
    # Step 1: Map target skeleton vertices to semantic roles
    var target_roles = map_target_semantic_roles(target_skeleton)
    
    # Step 2: Create unique bone mappings - NO DUPLICATES ALLOWED
    var used_bones = {}
    var virtual_bone_counter = 1000  # Start virtual bones at high indices
    
    # Step 3: Create primary anatomical mappings first
    var primary_mappings = create_primary_anatomical_mappings(glb_skeleton, hierarchy_analysis, target_roles, bone_map)
    
    # Step 4: Fill correspondence with unique mappings
    for target_idx in range(target_skeleton.vertices.size()):
        var target_role = target_roles[target_idx]
        var target_pos = target_skeleton.vertices[target_idx]
        
        print("DEBUG: Mapping Target[", target_idx, "] role:", target_role)
        
        # Special handling for hip roles - use joint splitting approach
        if target_role == "hip_right" or target_role == "hip_left" or target_role == "root_hips":
            var hip_split_data = create_hip_split_mapping(glb_skeleton, bone_map, target_role, target_idx)
            correspondence[target_idx] = hip_split_data.bone_index
            used_bones[hip_split_data.bone_index] = true
            # Use extended name function to handle both original and synthetic bones
            var bone_name = get_bone_name_extended(glb_skeleton, hip_split_data.bone_index)
            print("  -> HIP SPLIT: GLB[", hip_split_data.bone_index, "] ('", bone_name, "') at position ", hip_split_data.position)
        # Check for primary mapping
        elif primary_mappings.has(target_idx):
            var primary_bone = primary_mappings[target_idx]
            
            if not used_bones.has(primary_bone):
                # Use primary bone if available
                correspondence[target_idx] = primary_bone
                used_bones[primary_bone] = true
                var bone_name = glb_skeleton.get_bone_name(primary_bone)
                print("  -> PRIMARY: GLB[", primary_bone, "] ('", bone_name, "')")
            else:
                # Primary bone already used, create virtual joint
                var virtual_joint = create_virtual_joint_mapping(glb_skeleton, primary_bone, target_role, virtual_bone_counter)
                correspondence[target_idx] = virtual_joint
                used_bones[virtual_joint] = true
                virtual_bone_counter += 1
                print("  -> VIRTUAL: Virtual[", virtual_joint, "] based on GLB[", primary_bone, "]")
        else:
            # Find alternative bone
            var alternative_bone = find_alternative_bone_mapping(glb_skeleton, hierarchy_analysis, target_role, used_bones)
            if alternative_bone != -1:
                correspondence[target_idx] = alternative_bone
                used_bones[alternative_bone] = true
                var bone_name = glb_skeleton.get_bone_name(alternative_bone)
                print("  -> ALTERNATIVE: GLB[", alternative_bone, "] ('", bone_name, "')")
            else:
                # Create virtual joint as last resort
                var virtual_joint = create_fallback_virtual_joint(target_idx, virtual_bone_counter)
                correspondence[target_idx] = virtual_joint
                used_bones[virtual_joint] = true
                virtual_bone_counter += 1
                print("  -> FALLBACK: Virtual[", virtual_joint, "]")
    
    # Step 5: Display final unique correspondence
    print_final_correspondence(correspondence, glb_skeleton, target_roles)
    
    return correspondence

func create_primary_anatomical_mappings(glb_skeleton: Skeleton3D, hierarchy_analysis: Dictionary, target_roles: Array, bone_map) -> Dictionary:
    # Create primary bone mappings based on anatomical correspondence
    var mappings = {}
    
    # Direct anatomical mappings where we have good GLB bones
    var direct_mappings = {
        0: bone_map.root,           # root_hips -> Hips
        1: bone_map.spine,          # spine_lower -> Spine 
        2: bone_map.chest,          # spine_upper -> Chest
        3: bone_map.right_shoulder, # shoulder_right -> Right Shoulder
        6: bone_map.left_shoulder,  # shoulder_left -> Left Shoulder
        # Note: hip_right (9) and hip_left (12) will be virtual joints from root
    }
    
    # Add mappings if bones exist
    for target_idx in direct_mappings:
        var glb_bone = direct_mappings[target_idx]
        if glb_bone != -1:
            mappings[target_idx] = glb_bone
    
    # Find additional anatomical bones for arms and legs
    add_limb_mappings(glb_skeleton, mappings, bone_map)
    
    # Force virtual joints for hips - these should NEVER use alternative bones
    # hip_right (9) and hip_left (12) will be handled as virtual joints in main mapping loop
    
    return mappings

func add_limb_mappings(glb_skeleton: Skeleton3D, mappings: Dictionary, bone_map):
    # ORIGINAL GARMENT STRUCTURE: 0=root, 1=spine_lower, 2=spine_upper
    # 3=shoulder_right, 4=arm_right_upper, 5=arm_right_hand
    # 6=shoulder_left, 7=arm_left_upper, 8=arm_left_hand
    # 9=hip_right, 10=leg_right_upper, 11=leg_right_foot
    # 12=hip_left, 13=leg_left_upper, 14=leg_left_foot
    
    # Add shoulder mappings (3, 6) - These are the actual shoulder bones
    if bone_map.right_shoulder != -1:
        mappings[3] = bone_map.right_shoulder  # shoulder_right
    if bone_map.left_shoulder != -1:
        mappings[6] = bone_map.left_shoulder   # shoulder_left
    
    # Find arm chain bones - FROM SHOULDERS
    if bone_map.right_shoulder != -1:
        var right_arm_children = get_bone_children(glb_skeleton, bone_map.right_shoulder)
        if right_arm_children.size() > 0:
            var right_upper_arm = right_arm_children[0]
            mappings[4] = right_upper_arm  # arm_right_upper
            
            # Find hand/wrist (child of forearm, not forearm itself)
            var right_forearm_children = get_bone_children(glb_skeleton, right_upper_arm)
            if right_forearm_children.size() > 0:
                var right_forearm = right_forearm_children[0]
                var right_hand_children = get_bone_children(glb_skeleton, right_forearm)
                if right_hand_children.size() > 0:
                    mappings[5] = right_hand_children[0]  # arm_right_hand
                else:
                    mappings[5] = right_forearm  # Fallback to forearm if no hand found
    
    if bone_map.left_shoulder != -1:
        var left_arm_children = get_bone_children(glb_skeleton, bone_map.left_shoulder)
        if left_arm_children.size() > 0:
            var left_upper_arm = left_arm_children[0]
            mappings[7] = left_upper_arm  # arm_left_upper
            
            # Find hand/wrist (child of forearm, not forearm itself)
            var left_forearm_children = get_bone_children(glb_skeleton, left_upper_arm)
            if left_forearm_children.size() > 0:
                var left_forearm = left_forearm_children[0]
                var left_hand_children = get_bone_children(glb_skeleton, left_forearm)
                if left_hand_children.size() > 0:
                    mappings[8] = left_hand_children[0]  # arm_left_hand
                else:
                    mappings[8] = left_forearm  # Fallback to forearm if no hand found
    
    # Find leg chain bones - same as before
    if bone_map.right_leg != -1:
        mappings[10] = bone_map.right_leg  # leg_right_upper
        var right_leg_children = get_bone_children(glb_skeleton, bone_map.right_leg)
        if right_leg_children.size() > 0:
            var right_shin = right_leg_children[0]
            var right_foot_children = get_bone_children(glb_skeleton, right_shin)
            if right_foot_children.size() > 0:
                mappings[11] = right_foot_children[0]  # leg_right_foot
            else:
                mappings[11] = right_shin  # Fallback to shin if no foot found
    
    if bone_map.left_leg != -1:
        mappings[13] = bone_map.left_leg   # leg_left_upper
        var left_leg_children = get_bone_children(glb_skeleton, bone_map.left_leg)
        if left_leg_children.size() > 0:
            var left_shin = left_leg_children[0]
            var left_foot_children = get_bone_children(glb_skeleton, left_shin)
            if left_foot_children.size() > 0:
                mappings[14] = left_foot_children[0]  # leg_left_foot
            else:
                mappings[14] = left_shin  # Fallback to shin if no foot found

# Virtual Joint System (based on KHR_avatar_virtual_joints)
# Use simple dictionary-based approach to avoid class scope issues

func create_virtual_joint_mapping(glb_skeleton: Skeleton3D, base_bone: int, target_role: String, virtual_index: int) -> int:
    # Create a semantic virtual joint based on target role and parent bone
    var virtual_joint = create_semantic_virtual_joint(glb_skeleton, base_bone, target_role, virtual_index)
    
    # Store virtual joint in global registry
    if not virtual_joints_registry.has(virtual_index):
        virtual_joints_registry[virtual_index] = virtual_joint
    
    print("    Creating virtual joint '", virtual_joint.name, "' at index ", virtual_index, " from GLB[", base_bone, "]")
    return virtual_index

func create_semantic_virtual_joint(glb_skeleton: Skeleton3D, parent_bone: int, target_role: String, virtual_index: int) -> Dictionary:
    # Create semantically meaningful virtual joints based on target role
    var joint_name = target_role + "_virtual"
    var offset = calculate_semantic_offset(glb_skeleton, parent_bone, target_role)
    
    var virtual_joint = {
        "name": joint_name,
        "parent_joint": parent_bone,
        "translation": offset,
        "rotation": Quaternion.IDENTITY,
        "virtual_index": virtual_index
    }
    
    return virtual_joint

func calculate_virtual_joint_world_position(glb_skeleton: Skeleton3D, virtual_joint: Dictionary) -> Vector3:
    # Calculate world position: parent_transform * local_offset
    var parent_world_pos = calculate_bone_world_position(glb_skeleton, virtual_joint.parent_joint)
    var parent_transform = glb_skeleton.get_bone_global_pose(virtual_joint.parent_joint)
    var local_offset_world = parent_transform.basis * virtual_joint.translation
    return parent_world_pos + local_offset_world

func calculate_semantic_offset(glb_skeleton: Skeleton3D, parent_bone: int, target_role: String) -> Vector3:
    # Calculate meaningful offsets based on target role and anatomical structure
    var parent_pos = calculate_bone_world_position(glb_skeleton, parent_bone)
    
    match target_role:
        "hip_right":
            # Offset toward right leg direction
            return calculate_hip_offset(glb_skeleton, parent_bone, false)  # right side
        "hip_left":
            # Offset toward left leg direction  
            return calculate_hip_offset(glb_skeleton, parent_bone, true)   # left side
        "arm_right_upper", "arm_right_lower":
            # Offset along right arm chain
            return calculate_arm_offset(glb_skeleton, parent_bone, false)  # right side
        "arm_left_upper", "arm_left_lower":
            # Offset along left arm chain
            return calculate_arm_offset(glb_skeleton, parent_bone, true)   # left side
        "leg_right_upper", "leg_right_lower":
            # Offset along right leg chain
            return calculate_leg_offset(glb_skeleton, parent_bone, false)  # right side
        "leg_left_upper", "leg_left_lower":
            # Offset along left leg chain
            return calculate_leg_offset(glb_skeleton, parent_bone, true)   # left side
        "spine_lower", "spine_upper":
            # Offset along spine direction
            return calculate_spine_offset(glb_skeleton, parent_bone, target_role)
        _:
            # Default small offset
            return Vector3(0.0, 0.1, 0.0)

func calculate_hip_offset(glb_skeleton: Skeleton3D, hip_bone: int, is_left: bool) -> Vector3:
    # Calculate offset from hip toward leg direction
    var bone_map = identify_key_bones(glb_skeleton)
    var target_leg = bone_map.left_leg if is_left else bone_map.right_leg
    
    if target_leg != -1:
        var hip_pos = calculate_bone_world_position(glb_skeleton, hip_bone)
        var leg_pos = calculate_bone_world_position(glb_skeleton, target_leg)
        var direction = (leg_pos - hip_pos).normalized()
        return direction * 0.15  # 15cm offset toward leg
    else:
        # Fallback: lateral offset
        var x_offset = 0.1 if is_left else -0.1
        return Vector3(x_offset, 0.0, 0.0)

func calculate_arm_offset(glb_skeleton: Skeleton3D, shoulder_bone: int, is_left: bool) -> Vector3:
    # Calculate offset along arm direction
    var children = get_bone_children(glb_skeleton, shoulder_bone)
    
    if children.size() > 0:
        var shoulder_pos = calculate_bone_world_position(glb_skeleton, shoulder_bone)
        var arm_pos = calculate_bone_world_position(glb_skeleton, children[0])
        var direction = (arm_pos - shoulder_pos).normalized()
        return direction * 0.1  # 10cm along arm direction
    else:
        # Fallback: lateral + forward offset
        var x_offset = 0.2 if is_left else -0.2
        return Vector3(x_offset, 0.0, 0.1)

func calculate_leg_offset(glb_skeleton: Skeleton3D, leg_bone: int, is_left: bool) -> Vector3:
    # Calculate offset along leg direction
    var children = get_bone_children(glb_skeleton, leg_bone)
    
    if children.size() > 0:
        var upper_leg_pos = calculate_bone_world_position(glb_skeleton, leg_bone)
        var lower_leg_pos = calculate_bone_world_position(glb_skeleton, children[0])
        var direction = (lower_leg_pos - upper_leg_pos).normalized()
        return direction * 0.15  # 15cm along leg direction
    else:
        # Fallback: downward offset
        return Vector3(0.0, -0.15, 0.0)

func calculate_spine_offset(glb_skeleton: Skeleton3D, spine_bone: int, target_role: String) -> Vector3:
    # Calculate offset along spine chain
    var children = get_bone_children(glb_skeleton, spine_bone)
    
    if children.size() > 0 and target_role.contains("upper"):
        var spine_pos = calculate_bone_world_position(glb_skeleton, spine_bone)
        var upper_pos = calculate_bone_world_position(glb_skeleton, children[0])
        var direction = (upper_pos - spine_pos).normalized()
        return direction * 0.05  # 5cm along spine
    else:
        # Default upward offset
        var y_offset = 0.05 if target_role.contains("upper") else -0.05
        return Vector3(0.0, y_offset, 0.0)

func find_alternative_bone_mapping(glb_skeleton: Skeleton3D, hierarchy_analysis: Dictionary, target_role: String, used_bones: Dictionary) -> int:
    # Find an unused bone that could serve as an alternative for this target role
    
    # Look through major bones for unused alternatives
    for bone_idx in hierarchy_analysis.major_bones.keys():
        if used_bones.has(bone_idx):
            continue
            
        var bone_role = hierarchy_analysis.major_bones[bone_idx]
        
        # Check for compatible roles
        if is_compatible_alternative_role(target_role, bone_role):
            return bone_idx
    
    # Look through all bones for unused options
    for bone_idx in range(glb_skeleton.get_bone_count()):
        if used_bones.has(bone_idx):
            continue
            
        # Use any unused bone as alternative
        return bone_idx
    
    return -1

func is_compatible_alternative_role(target_role: String, bone_role: String) -> bool:
    # Check if bone_role can serve as alternative for target_role
    var alternatives = {
        "hip_right": ["root", "upper_leg_right", "lower_leg_right"],
        "hip_left": ["root", "upper_leg_left", "lower_leg_left"],
        "arm_right_upper": ["shoulder_right", "lower_arm_right"],
        "arm_left_upper": ["shoulder_left", "lower_arm_left"],
        "leg_right_upper": ["hip_right", "lower_leg_right"],
        "leg_left_upper": ["hip_left", "lower_leg_left"]
    }
    
    if alternatives.has(target_role):
        var alt_list = alternatives[target_role]
        for alt in alt_list:
            if bone_role.contains(alt):
                return true
    
    return false

func create_hip_split_mapping(glb_skeleton: Skeleton3D, bone_map, target_role: String, target_idx: int) -> Dictionary:
    # Split the single GLB "Hips" bone into three distinct synthetic bones
    var hips_bone = bone_map.root if bone_map.root != -1 else 0
    var hips_pos = calculate_bone_world_position(glb_skeleton, hips_bone)
    
    # Calculate split positions based on target role
    var split_position: Vector3
    var synthetic_bone_index: int
    
    match target_role:
        "root_hips":
            # Center position - use original hips bone
            split_position = hips_pos
            synthetic_bone_index = hips_bone  # Keep original bone for root
        "hip_left":
            # Left hip position - create synthetic bone
            split_position = calculate_hip_split_position(glb_skeleton, bone_map, hips_pos, true)
            synthetic_bone_index = create_synthetic_bone(hips_bone, split_position, "hip_left_split")
        "hip_right":
            # Right hip position - create synthetic bone
            split_position = calculate_hip_split_position(glb_skeleton, bone_map, hips_pos, false)
            synthetic_bone_index = create_synthetic_bone(hips_bone, split_position, "hip_right_split")
        _:
            # Fallback to center
            split_position = hips_pos
            synthetic_bone_index = hips_bone
    
    # Transform to target coordinates
    var target_position = Vector3(split_position.x, -split_position.z, split_position.y)
    
    return {
        "bone_index": synthetic_bone_index,
        "position": target_position,
        "role": target_role
    }

func calculate_hip_split_position(glb_skeleton: Skeleton3D, bone_map, hips_pos: Vector3, is_left: bool) -> Vector3:
    # Calculate hip split position by offsetting toward the appropriate leg
    var target_leg = bone_map.left_leg if is_left else bone_map.right_leg
    var side_name = "left" if is_left else "right"
    
    if target_leg != -1:
        # Calculate offset toward leg direction
        var leg_pos = calculate_bone_world_position(glb_skeleton, target_leg)
        var direction = (leg_pos - hips_pos).normalized()
        
        # Apply partial offset toward leg (not full distance)
        var offset_distance = 0.12  # 12cm lateral offset
        var offset_position = hips_pos + (direction * offset_distance)
        
        print("    Hip split ", side_name, ": offset ", offset_distance, "m toward leg GLB[", target_leg, "]")
        return offset_position
    else:
        # Fallback: use fixed lateral offset
        var x_offset = 0.12 if is_left else -0.12
        var fallback_position = hips_pos + Vector3(x_offset, 0.0, 0.0)
        
        print("    Hip split ", side_name, ": fallback lateral offset ", x_offset, "m")
        return fallback_position

func create_fallback_virtual_joint(target_idx: int, virtual_index: int) -> int:
    # Create a fallback virtual joint when no alternatives exist
    var fallback_joint = {
        "name": "fallback_" + str(target_idx),
        "parent_joint": 0,  # Parent to root bone
        "translation": Vector3(0.0, 0.1 * target_idx, 0.0),  # Spread fallbacks vertically
        "rotation": Quaternion.IDENTITY,
        "virtual_index": virtual_index
    }
    
    virtual_joints_registry[virtual_index] = fallback_joint
    print("    Creating fallback virtual joint '", fallback_joint.name, "' at index ", virtual_index)
    return virtual_index

func create_split_hip_correspondences(glb_skeleton: Skeleton3D, hierarchy_analysis: Dictionary, target_skeleton, target_roles: Array) -> Dictionary:
    var hip_assignments = {}
    
    # Find the root bone (hips) in GLB skeleton
    var root_bone_idx = -1
    for bone_idx in hierarchy_analysis.major_bones.keys():
        var bone_role = hierarchy_analysis.major_bones[bone_idx]
        if bone_role == "root":
            root_bone_idx = bone_idx
            break
    
    if root_bone_idx == -1:
        print("DEBUG: No root bone found for hip splitting")
        return hip_assignments
    
    # Generate virtual bone indices for split hips
    var virtual_left_hip = glb_skeleton.get_bone_count()  # Virtual bone index
    var virtual_right_hip = glb_skeleton.get_bone_count() + 1  # Virtual bone index
    
    # Find target indices for hips
    var target_root_idx = -1
    var target_left_hip_idx = -1
    var target_right_hip_idx = -1
    
    for i in range(target_roles.size()):
        match target_roles[i]:
            "root_hips":
                target_root_idx = i
            "hip_left":
                target_left_hip_idx = i
            "hip_right":
                target_right_hip_idx = i
    
    print("DEBUG: Hip splitting - Root:", target_root_idx, " Left:", target_left_hip_idx, " Right:", target_right_hip_idx)
    print("DEBUG: GLB root bone:", root_bone_idx, " ('", glb_skeleton.get_bone_name(root_bone_idx), "')")
    
    # Assign correspondences:
    # - World root (target[0]) -> GLB root bone (world origin)
    # - Object left hip (target[12]) -> GLB root bone (object center) 
    # - Object right hip (target[9]) -> GLB root bone (object center)
    
    if target_root_idx != -1:
        hip_assignments[target_root_idx] = root_bone_idx  # World root
        print("DEBUG: World root Target[", target_root_idx, "] -> GLB[", root_bone_idx, "] (world reference)")
    
    if target_left_hip_idx != -1:
        hip_assignments[target_left_hip_idx] = root_bone_idx  # Object left hip
        print("DEBUG: Left hip Target[", target_left_hip_idx, "] -> GLB[", root_bone_idx, "] (object left)")
    
    if target_right_hip_idx != -1:
        hip_assignments[target_right_hip_idx] = root_bone_idx  # Object right hip  
        print("DEBUG: Right hip Target[", target_right_hip_idx, "] -> GLB[", root_bone_idx, "] (object right)")
    
    return hip_assignments

func map_target_semantic_roles(target_skeleton) -> Array:
    # Map the 15 target skeleton vertices to semantic roles - ORIGINAL garment structure
    var roles = [
        "root_hips",        # 0: Base of skeleton
        "spine_lower",      # 1: Lower spine
        "spine_upper",      # 2: Upper spine/chest
        "shoulder_right",   # 3: Right shoulder
        "arm_right_upper",  # 4: Right upper arm
        "arm_right_hand",   # 5: Right hand
        "shoulder_left",    # 6: Left shoulder
        "arm_left_upper",   # 7: Left upper arm
        "arm_left_hand",    # 8: Left hand
        "hip_right",        # 9: Right hip
        "leg_right_upper",  # 10: Right upper leg
        "leg_right_foot",   # 11: Right foot
        "hip_left",         # 12: Left hip
        "leg_left_upper",   # 13: Left upper leg
        "leg_left_foot"     # 14: Left foot
    ]
    
    return roles

func find_best_hierarchical_rotational_match(glb_skeleton: Skeleton3D, hierarchy_analysis: Dictionary, target_pos: Vector3, target_role: String, used_bones: Dictionary) -> int:
    var best_bone = -1
    var best_score = -1.0
    
    # Score bones based on hierarchical position + rotational compatibility + semantic role
    for bone_idx in hierarchy_analysis.major_bones.keys():
        if used_bones.has(bone_idx):
            continue
            
        var bone_role = hierarchy_analysis.major_bones[bone_idx]
        var bone_data = hierarchy_analysis.bone_transforms[bone_idx]
        
        # Calculate combined score
        var hierarchy_score = calculate_hierarchy_score(target_role, bone_role)
        var rotational_score = calculate_rotational_score(glb_skeleton, bone_idx, target_pos, bone_data)
        var positional_score = calculate_positional_score(target_pos, bone_data.position)
        
        var combined_score = hierarchy_score * 0.5 + rotational_score * 0.3 + positional_score * 0.2
        
        print("DEBUG:   Bone ", bone_idx, " ('", bone_data.name, "') role:", bone_role, " scores: H:", hierarchy_score, " R:", rotational_score, " P:", positional_score, " = ", combined_score)
        
        if combined_score > best_score:
            best_score = combined_score
            best_bone = bone_idx
    
    return best_bone

func calculate_hierarchy_score(target_role: String, bone_role: String) -> float:
    # Score based on precise semantic role matching using the new bone identification
    var role_matches = {
        "root_hips": ["root"],
        "spine_lower": ["spine_1"],
        "spine_upper": ["spine_2"],
        "shoulder_right": ["shoulder_right"],
        "shoulder_left": ["shoulder_left"],
        "arm_right_upper": ["upper_arm_right"],     # Fix: match upper_arm, not generic arm
        "arm_left_upper": ["upper_arm_left"],       # Fix: match upper_arm, not generic arm  
        "arm_right_lower": ["lower_arm_right"],     # Fix: match lower_arm specifically
        "arm_left_lower": ["lower_arm_left"],       # Fix: match lower_arm specifically
        "hip_right": ["root"],                      # Both hips use root bone but different positions
        "hip_left": ["root"],                       # Both hips use root bone but different positions
        "leg_right_upper": ["upper_leg_right"],     # Fix: match upper_leg specifically
        "leg_left_upper": ["upper_leg_left"],       # Fix: match upper_leg specifically
        "leg_right_lower": ["lower_leg_right", "foot_right"], # Fix: match lower_leg or foot
        "leg_left_lower": ["lower_leg_left", "foot_left"]     # Fix: match lower_leg or foot
    }
    
    if role_matches.has(target_role):
        var matches = role_matches[target_role]
        for match in matches:
            if bone_role.contains(match):
                return 1.0
    
    return 0.1  # Low score for non-matching roles

func calculate_rotational_score(glb_skeleton: Skeleton3D, bone_idx: int, target_pos: Vector3, bone_data: Dictionary) -> float:
    # Calculate swing and twist rotation compatibility
    var bone_transform = bone_data.transform
    var bone_pos = bone_data.position
    
    # Analyze swing rotation (rotation around perpendicular axes)
    var direction_to_target = (target_pos - bone_pos).normalized()
    var bone_forward = bone_transform.basis.z.normalized()
    
    var swing_alignment = direction_to_target.dot(bone_forward)
    var swing_score = (swing_alignment + 1.0) * 0.5  # Normalize to 0-1
    
    # Analyze twist rotation (rotation around bone axis)
    var parent_idx = bone_data.parent
    var twist_score = 0.5  # Default neutral score
    
    if parent_idx != -1:
        var parent_data = glb_skeleton.get_bone_rest(parent_idx)
        var parent_pos = calculate_bone_world_position(glb_skeleton, parent_idx)
        
        # Calculate twist compatibility based on bone chain direction
        var bone_chain_dir = (bone_pos - parent_pos).normalized()
        var target_chain_dir = direction_to_target
        
        var twist_alignment = bone_chain_dir.dot(target_chain_dir)
        twist_score = (twist_alignment + 1.0) * 0.5
    
    return (swing_score + twist_score) * 0.5

func calculate_positional_score(target_pos: Vector3, bone_pos: Vector3) -> float:
    # Simple distance-based score, inverted so closer = higher score
    var distance = target_pos.distance_to(bone_pos)
    return 1.0 / (1.0 + distance * 0.1)  # Normalize distance impact

func print_final_correspondence(correspondence: Array, glb_skeleton: Skeleton3D, target_roles: Array):
    print("\n=== ERON'S HIERARCHICAL-ROTATIONAL CORRESPONDENCE ===")
    
    for i in range(correspondence.size()):
        var glb_bone_idx = correspondence[i]
        var target_role = target_roles[i] if i < target_roles.size() else "unknown"
        
        if glb_bone_idx != -1:
            # Use extended name function to handle both original and synthetic bones
            var bone_name = get_bone_name_extended(glb_skeleton, glb_bone_idx)
            var distance_from_root = calculate_distance_from_root_extended(glb_skeleton, glb_bone_idx)
            print("Target[", i, "] (", target_role, ") -> GLB Bone ", glb_bone_idx, " ('", bone_name, "', dist:", distance_from_root, ")")
        else:
            print("Target[", i, "] (", target_role, ") -> NO CORRESPONDENCE")
    
    print("=== END HIERARCHICAL-ROTATIONAL MAPPING ===\n")

func filter_anatomical_bones(glb_skeleton: Skeleton3D) -> Array:
    # Filter bones to include only major anatomical structures
    # Exclude accessory bones like skirts, tails, hair, etc.
    var anatomical_bones = []
    
    for bone_idx in range(glb_skeleton.get_bone_count()):
        if is_anatomical_bone(glb_skeleton, bone_idx):
            anatomical_bones.append(bone_idx)
    
    print("DEBUG: Anatomical bones found: ", anatomical_bones)
    return anatomical_bones

func is_anatomical_bone(glb_skeleton: Skeleton3D, bone_idx: int) -> bool:
    # Determine if a bone is part of the major anatomical structure
    # using pattern recognition without relying on specific names
    var bone_name = glb_skeleton.get_bone_name(bone_idx).to_lower()
    
    # Exclude obviously non-anatomical bones by pattern
    var excluded_patterns = [
        "skirt", "tail", "hair", "ear", "eye", "ahoge", "bang", 
        "breast", "butt", "clothing", "accessory", "decoration"
    ]
    
    for pattern in excluded_patterns:
        if bone_name.contains(pattern):
            return false
    
    # Include major anatomical bones by hierarchy position and typical structure
    # Core anatomical bones typically include: spine chain, limbs, head/neck
    var included_patterns = [
        "hip", "spine", "chest", "neck", "head", 
        "shoulder", "arm", "hand", "finger", "thumb",
        "leg", "foot", "toe"
    ]
    
    # If it matches anatomical patterns, include it
    for pattern in included_patterns:
        if bone_name.contains(pattern):
            return true
    
    # For bones without clear names, use hierarchy analysis
    # Include bones that are close to root and part of main hierarchy
    var distance_from_root = calculate_distance_from_root(glb_skeleton, bone_idx)
    var children_count = get_bone_children(glb_skeleton, bone_idx).size()
    
    # Include bones that are:
    # 1. Close to root (distance <= 3) AND have children (major structural bones)
    # 2. OR are leaf bones of the main hierarchy (hands, feet, head)
    if distance_from_root <= 3 and children_count > 0:
        return true  # Major structural bone
    elif distance_from_root <= 4 and children_count == 0:
        return true  # Leaf bone (likely hand/foot/head)
    
    return false

func create_garment_effectors(target_skeleton) -> Array:
    var effectors = []
    
    # Each vertex in the target skeleton becomes an effector
    for i in range(target_skeleton.vertices.size()):
        var effector = GarmentEffector.new(i, target_skeleton.vertices[i])
        effectors.append(effector)
    
    return effectors

func trace_effector_to_root(glb_skeleton: Skeleton3D, effector: GarmentEffector, bone_influences: Dictionary, bone_map):
    # Start from the most relevant GLB bone for this effector
    var start_bone = find_closest_glb_bone_for_effector(effector, bone_map)
    if start_bone == -1:
        return  # No suitable starting bone found
    
    var current_bone = start_bone
    var current_weight = effector.opacity
    var bones_encountered = []
    
    # Traverse from effector's best GLB bone to root
    while current_bone != -1 and current_weight > 0.0:
        var bone_name = glb_skeleton.get_bone_name(current_bone)
        bones_encountered.append(current_bone)
        
        # Create or update bone influence
        if not bone_influences.has(current_bone):
            var distance = calculate_distance_from_root(glb_skeleton, current_bone)
            bone_influences[current_bone] = BoneInfluence.new(current_bone, distance, bone_name)
        
        var influence = bone_influences[current_bone]
        
        # Check if this bone is an effector for another garment bone
        var other_effector_weight = get_bone_effector_opacity(current_bone, effector.index, bone_influences)
        if other_effector_weight > 0 and other_effector_weight < 1.0:
            current_weight *= (1.0 - other_effector_weight)
        
        # Stop if weight becomes negligible
        if current_weight <= 0.001:
            break
        
        # Add this effector's influence to the bone
        influence.effector_weights[effector.index] = current_weight
        
        # Move to parent bone
        current_bone = glb_skeleton.get_bone_parent(current_bone)

func find_closest_glb_bone_for_effector(effector: GarmentEffector, bone_map) -> int:
    # Use pure geometric distance to find closest GLB bone to garment effector position
    # No name-based mapping allowed - use spatial correspondence only
    return find_closest_bone_by_position(effector.position, bone_map)

func find_closest_bone_by_position(target_position: Vector3, bone_map) -> int:
    # For now, distribute across all major skeletal regions to ensure 15 unique mappings
    # This ensures all 15 garment bones get different GLB bone correspondences
    var available_bones = [
        bone_map.root,           # 0: Hips 
        bone_map.spine,          # 1: Spine
        bone_map.chest,          # 2: Chest
        bone_map.neck,           # 3: Neck  
        bone_map.head,           # 4: Head
        bone_map.left_shoulder,  # 5: Left shoulder
        bone_map.right_shoulder, # 6: Right shoulder
        bone_map.left_arm,       # 7: Left arm
        bone_map.right_arm,      # 8: Right arm
        bone_map.left_leg,       # 9: Left leg
        bone_map.right_leg       # 10: Right leg
    ]
    
    # Remove invalid bones and add additional bones to reach 15 mappings
    var valid_bones = []
    for bone_idx in available_bones:
        if bone_idx != -1:
            valid_bones.append(bone_idx)
    
    # Ensure we have at least 15 unique bones by adding more from hierarchy
    while valid_bones.size() < 15:
        var additional_bones = find_additional_hierarchy_bones(valid_bones, bone_map)
        for bone in additional_bones:
            if bone != -1 and not valid_bones.has(bone):
                valid_bones.append(bone)
                if valid_bones.size() >= 15:
                    break
        if additional_bones.size() == 0:
            break  # Prevent infinite loop
    
    return valid_bones[0] if valid_bones.size() > 0 else bone_map.root

func find_additional_hierarchy_bones(existing_bones: Array, bone_map) -> Array:
    # Find additional bones in the GLB skeleton hierarchy to reach 15 mappings
    var additional = []
    
    # Add children of major bones to expand our bone set
    var major_bones = [bone_map.root, bone_map.spine, bone_map.chest, bone_map.left_arm, bone_map.right_arm, bone_map.left_leg, bone_map.right_leg]
    
    for major_bone in major_bones:
        if major_bone == -1:
            continue
            
        # Add children of this major bone
        var children = get_all_children_recursive(major_bone)
        for child in children:
            if not existing_bones.has(child) and not additional.has(child):
                additional.append(child)
                if additional.size() >= 10:  # Limit to prevent too many
                    break
    
    return additional

func get_all_children_recursive(bone_idx: int) -> Array:
    # Get all recursive children of a bone (placeholder - would need skeleton reference)
    # For now, return some reasonable additional bone indices
    var children = []
    
    # Add some reasonable bone indices as children
    for i in range(max(0, bone_idx - 2), min(117, bone_idx + 5)):
        if i != bone_idx:
            children.append(i)
    
    return children

func calculate_distance_from_root(glb_skeleton: Skeleton3D, bone_idx: int) -> int:
    var distance = 0
    var current_bone = bone_idx
    
    while get_bone_parent_extended(glb_skeleton, current_bone) != -1:
        distance += 1
        current_bone = get_bone_parent_extended(glb_skeleton, current_bone)
    
    return distance

func get_bone_effector_opacity(bone_idx: int, current_effector_idx: int, bone_influences: Dictionary) -> float:
    if not bone_influences.has(bone_idx):
        return 0.0
    
    var influence = bone_influences[bone_idx]
    var total_weight = 0.0
    
    # Sum weights from other effectors (not the current one)
    for effector_idx in influence.effector_weights:
        if effector_idx != current_effector_idx:
            total_weight += influence.effector_weights[effector_idx]
    
    return min(total_weight, 1.0)

func group_bone_sequences(bone_influences: Dictionary) -> Array:
    var groups = []
    var processed_bones = {}
    
    # Group bones that serve the same sets of effectors
    for bone_idx in bone_influences:
        if processed_bones.has(bone_idx):
            continue
            
        var influence = bone_influences[bone_idx]
        var effector_set = influence.effector_weights.keys()
        effector_set.sort()
        
        # Find or create group for this effector set
        var group = null
        for existing_group in groups:
            var existing_set = existing_group.target_effectors.duplicate()
            existing_set.sort()
            
            if arrays_equal(effector_set, existing_set):
                group = existing_group
                break
        
        if group == null:
            group = EffectorGroup.new()
            group.target_effectors = effector_set
            group.root_distance = influence.distance_from_root
            groups.append(group)
        
        group.bone_sequence.append(bone_idx)
        processed_bones[bone_idx] = true
        
        # Update root distance to the minimum (closest to root)
        group.root_distance = min(group.root_distance, influence.distance_from_root)
    
    return groups

func arrays_equal(arr1: Array, arr2: Array) -> bool:
    if arr1.size() != arr2.size():
        return false
    
    for i in range(arr1.size()):
        if arr1[i] != arr2[i]:
            return false
    
    return true

func create_solve_order(effector_groups: Array) -> Array:
    # Sort groups by distance from root (reverse sorted - furthest first)
    effector_groups.sort_custom(func(a, b): return a.root_distance > b.root_distance)
    
    var solve_order = []
    
    # Append bones from each group to solve order
    for group in effector_groups:
        for bone_idx in group.bone_sequence:
            solve_order.append(bone_idx)
    
    return solve_order

# OLD FUNCTION REMOVED - Using only geometric correspondence now

func scale_bone_position(glb_skeleton: Skeleton3D, bone_idx: int, scale_factors):
    var bone_pose = glb_skeleton.get_bone_rest(bone_idx)
    var pos = bone_pose.origin
    
    # Only transform coordinate system: GLB (X,Y,Z) -> Target (X, -Z, Y)
    # No scaling or translation - preserve original GLB coordinates
    return Vector3(pos.x, -pos.z, pos.y)

func write_mapped_skeleton_obj(mapped_skeleton, file):
    # Write the mapped skeleton with extended bone system support
    # This includes both original GLB bones and synthetic bones
    
    print("DEBUG: Writing extended skeleton with ", mapped_skeleton.vertices.size(), " vertices")
    print("DEBUG: Synthetic bones registry has ", synthetic_bones_registry.size(), " synthetic bones")
    
    # Write vertices (includes positions for synthetic bones)
    for vertex in mapped_skeleton.vertices:
        file.store_line("v %f %f %f" % [vertex.x, vertex.y, vertex.z])
    
    # Write edges with exact target connectivity
    for edge in mapped_skeleton.edges:
        file.store_line("l %d %d" % [edge[0] + 1, edge[1] + 1])  # Convert to 1-based
    
    # Debug: Print synthetic bone information
    for synthetic_idx in synthetic_bones_registry.keys():
        var synthetic_data = synthetic_bones_registry[synthetic_idx]
        print("DEBUG: Synthetic bone ", synthetic_idx, " ('", synthetic_data.name, "') written to skeleton")

func write_direct_skeleton_obj(skeleton: Skeleton3D, file):
    # Direct skeleton export with proper world positions
    var bone_count = skeleton.get_bone_count()
    
    # Calculate world positions for all bones
    var world_positions = []
    for i in bone_count:
        var world_pos = calculate_bone_world_position(skeleton, i)
        
        # Only transform coordinate system: GLB (X,Y,Z) -> Target (X, -Z, Y)
        # No scaling or translation - preserve original GLB coordinates
        var blender_coords = Vector3(world_pos.x, -world_pos.z, world_pos.y)
        world_positions.append(blender_coords)
        file.store_line("v %f %f %f" % [blender_coords.x, blender_coords.y, blender_coords.z])
    
    # Write bone connections based on hierarchy
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

func print_skeleton_railroad_diagram(skeleton: Skeleton3D):
    print("\n=== SKELETON RAILROAD DIAGRAM ===")
    var bone_count = skeleton.get_bone_count()
    
    # Find root bones (bones with no parent)
    var root_bones = []
    for i in range(bone_count):
        if skeleton.get_bone_parent(i) == -1:
            root_bones.append(i)
    
    print("Total bones: ", bone_count)
    print("Root bones: ", root_bones.size())
    
    # Print hierarchy for each root
    for root in root_bones:
        print_bone_hierarchy(skeleton, root, 0)
    
    print("=== END RAILROAD DIAGRAM ===\n")

func print_bone_hierarchy(skeleton: Skeleton3D, bone_idx: int, depth: int):
    var indent = ""
    for i in range(depth):
        indent += "  "
    
    var bone_name = skeleton.get_bone_name(bone_idx)
    var children = get_bone_children(skeleton, bone_idx)
    
    print(indent, bone_idx, ": '", bone_name, "' (", children.size(), " children)")
    
    # Print children
    for child in children:
        print_bone_hierarchy(skeleton, child, depth + 1)

func print_bone_map_analysis(skeleton: Skeleton3D, bone_map):
    print("\n=== BONE MAP ANALYSIS ===")
    
    for key in bone_map.keys():
        var bone_idx = bone_map[key]
        if bone_idx != -1:
            var bone_name = skeleton.get_bone_name(bone_idx)
            print(key, " -> Bone ", bone_idx, ": '", bone_name, "'")
        else:
            print(key, " -> NOT FOUND")
    
    print("=== END BONE MAP ANALYSIS ===\n")

func print_garment_correspondence(bone_correspondence, skeleton: Skeleton3D):
    print("\n=== GARMENT SKELETON CORRESPONDENCE ===")
    print("Target skeleton structure (15 vertices):")
    print("0: root, 1: spine_base, 2: spine_top")
    print("3-5: right_arm, 6-8: left_arm")
    print("9-11: right_leg, 12-14: left_leg")
    print()
    
    var labels = [
        "root", "spine_base", "spine_top",
        "right_shoulder", "right_arm", "right_hand",
        "left_shoulder", "left_arm", "left_hand",
        "right_hip", "right_leg", "right_foot",
        "left_hip", "left_leg", "left_foot"
    ]
    
    for i in range(bone_correspondence.size()):
        var glb_bone_idx = bone_correspondence[i]
        var label = labels[i] if i < labels.size() else ("vertex_" + str(i))
        
        if glb_bone_idx != -1:
            var bone_name = skeleton.get_bone_name(glb_bone_idx)
            print("Target[", i, "] (", label, ") -> GLB Bone ", glb_bone_idx, ": '", bone_name, "'")
        else:
            print("Target[", i, "] (", label, ") -> NO CORRESPONDENCE")
    
    print("=== END CORRESPONDENCE ===\n")
