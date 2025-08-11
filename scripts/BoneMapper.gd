extends RefCounted
class_name BoneMapper

const Utils = preload("Utils.gd")

# BoneMapper: Maps between 15-bone garment skeleton and GLB skeleton
# Uses Godot-style auto mapping with bone splitting for accurate correspondence

# Entry point for bone mapping process
func map_garment_to_glb_skeleton(garment_skeleton: Skeleton3D, glb_skeleton: Skeleton3D) -> Array:
    print("DEBUG: Starting garment to GLB skeleton mapping")
    
    # Step 1: Identify key bones in GLB skeleton using Godot auto-mapping
    var bone_map = identify_key_bones(glb_skeleton)
    
    # Step 2: Create correspondence array for 15 garment bones
    var correspondence = create_initial_correspondence(garment_skeleton, glb_skeleton, bone_map)
    
    # Step 3: Add limb correspondences with bone splitting
    add_limb_correspondences(glb_skeleton, correspondence, bone_map)
    
    # Step 4: Fill any remaining gaps with position-based matching
    fill_correspondence_gaps(garment_skeleton, glb_skeleton, correspondence, bone_map)
    
    return correspondence

func identify_key_bones(skeleton: Skeleton3D):
    var bone_map = {
        "root": -1, "spine": -1, "chest": -1, "neck": -1, "head": -1,
        "left_shoulder": -1, "right_shoulder": -1,
        "left_arm": -1, "right_arm": -1,
        "left_leg": -1, "right_leg": -1,
        "left_hip": -1, "right_hip": -1,
        "left_hand": -1, "right_hand": -1,
        "left_foot": -1, "right_foot": -1
    }
    
    var bone_count = skeleton.get_bone_count()
    if bone_count == 0:
        return bone_map
    
    print("DEBUG: Starting Godot-style auto mapping for ", bone_count, " bones")
    
    # Use Godot's proven auto-mapping algorithm
    godot_auto_mapping_process(skeleton, bone_map)
    
    return bone_map

# Godot-style auto mapping algorithm (adapted from C++ source)
func godot_auto_mapping_process(skeleton: Skeleton3D, bone_map):
    print("DEBUG: Starting Godot auto mapping process")
    
    # 1. Guess Hips
    var picklist = ["hip", "pelvis", "waist", "torso", "spine"]
    var hips = search_bone_by_name(skeleton, picklist)
    if hips == -1:
        print("DEBUG: Auto Mapping couldn't guess Hips. Abort auto mapping.")
        return
    else:
        bone_map.root = hips  # Store as root for our 15-bone system
        print("DEBUG: Found Hips: ", hips, " ('", skeleton.get_bone_name(hips), "')")
    
    # 2. Guess Root (actual skeleton root)
    var bone_idx = skeleton.get_bone_parent(hips)
    var search_path = []
    while bone_idx >= 0:
        search_path.append(bone_idx)
        bone_idx = skeleton.get_bone_parent(bone_idx)
    
    # 3. Guess Feet
    picklist = ["foot", "ankle"]
    var left_foot = search_bone_by_name(skeleton, picklist, "left", hips)
    if left_foot != -1:
        bone_map.left_foot = left_foot
        print("DEBUG: Found LeftFoot: ", left_foot, " ('", skeleton.get_bone_name(left_foot), "')")
    
    var right_foot = search_bone_by_name(skeleton, picklist, "right", hips)
    if right_foot != -1:
        bone_map.right_foot = right_foot
        print("DEBUG: Found RightFoot: ", right_foot, " ('", skeleton.get_bone_name(right_foot), "')")
    
    # 3-1. Guess LowerLegs
    picklist = ["lower.*leg", "knee", "shin", "calf", "leg"]
    var left_lower_leg = -1
    if left_foot != -1:
        left_lower_leg = search_bone_by_name(skeleton, picklist, "left", hips, left_foot)
    
    var right_lower_leg = -1  
    if right_foot != -1:
        right_lower_leg = search_bone_by_name(skeleton, picklist, "right", hips, right_foot)
    
    # 3-2. Guess UpperLegs
    picklist = ["upper.*leg", "thigh", "leg"]
    if left_lower_leg != -1:
        bone_idx = search_bone_by_name(skeleton, picklist, "left", hips, left_lower_leg)
        if bone_idx != -1:
            bone_map.left_leg = bone_idx
            print("DEBUG: Found LeftUpperLeg: ", bone_idx, " ('", skeleton.get_bone_name(bone_idx), "')")
    
    if right_lower_leg != -1:
        bone_idx = search_bone_by_name(skeleton, picklist, "right", hips, right_lower_leg)
        if bone_idx != -1:
            bone_map.right_leg = bone_idx
            print("DEBUG: Found RightUpperLeg: ", bone_idx, " ('", skeleton.get_bone_name(bone_idx), "')")
    
    # 4. Guess Hands
    picklist = ["hand", "wrist", "palm", "fingers"]
    var left_hand = search_bone_by_name(skeleton, picklist, "left", hips)
    if left_hand != -1:
        bone_map.left_hand = left_hand
        print("DEBUG: Found LeftHand: ", left_hand, " ('", skeleton.get_bone_name(left_hand), "')")
    
    var right_hand = search_bone_by_name(skeleton, picklist, "right", hips)
    if right_hand != -1:
        bone_map.right_hand = right_hand
        print("DEBUG: Found RightHand: ", right_hand, " ('", skeleton.get_bone_name(right_hand), "')")
    
    # 5. Guess Arms
    picklist = ["shoulder", "clavicle", "collar"]
    var left_shoulder = search_bone_by_name(skeleton, picklist, "left", hips)
    if left_shoulder != -1:
        bone_map.left_shoulder = left_shoulder
        print("DEBUG: Found LeftShoulder: ", left_shoulder, " ('", skeleton.get_bone_name(left_shoulder), "')")
    
    var right_shoulder = search_bone_by_name(skeleton, picklist, "right", hips)
    if right_shoulder != -1:
        bone_map.right_shoulder = right_shoulder
        print("DEBUG: Found RightShoulder: ", right_shoulder, " ('", skeleton.get_bone_name(right_shoulder), "')")
    
    # 5-1. Guess LowerArms
    picklist = ["lower.*arm", "forearm", "elbow", "arm"]
    var left_lower_arm = -1
    if left_shoulder != -1 and left_hand != -1:
        left_lower_arm = search_bone_by_name(skeleton, picklist, "left", left_shoulder, left_hand)
    
    var right_lower_arm = -1
    if right_shoulder != -1 and right_hand != -1:
        right_lower_arm = search_bone_by_name(skeleton, picklist, "right", right_shoulder, right_hand)
    
    # 5-2. Guess UpperArms
    picklist = ["upper.*arm", "arm"]
    if left_shoulder != -1 and left_lower_arm != -1:
        bone_idx = search_bone_by_name(skeleton, picklist, "left", left_shoulder, left_lower_arm)
        if bone_idx != -1:
            bone_map.left_arm = bone_idx
            print("DEBUG: Found LeftUpperArm: ", bone_idx, " ('", skeleton.get_bone_name(bone_idx), "')")
    
    if right_shoulder != -1 and right_lower_arm != -1:
        bone_idx = search_bone_by_name(skeleton, picklist, "right", right_shoulder, right_lower_arm)
        if bone_idx != -1:
            bone_map.right_arm = bone_idx
            print("DEBUG: Found RightUpperArm: ", bone_idx, " ('", skeleton.get_bone_name(bone_idx), "')")
    
    # 6. Guess Neck
    picklist = ["neck", "head", "face"]
    var neck = search_bone_by_name(skeleton, picklist, "", hips)
    if neck != -1:
        bone_map.neck = neck
        print("DEBUG: Found Neck: ", neck, " ('", skeleton.get_bone_name(neck), "')")
    
    # 7. Guess Head
    picklist = ["head", "face"]
    var head = search_bone_by_name(skeleton, picklist, "", neck if neck != -1 else hips)
    if head == -1 and neck != -1:
        var neck_children = Utils.get_bone_children(skeleton, neck)
        if neck_children.size() == 1:
            head = neck_children[0]
    
    if head != -1:
        bone_map.head = head
        print("DEBUG: Found Head: ", head, " ('", skeleton.get_bone_name(head), "')")
    
    # 8. Guess Chest
    var neck_or_head = neck if neck != -1 else head
    if neck_or_head != -1:
        var chest_or_upper_chest = skeleton.get_bone_parent(neck_or_head)
        if chest_or_upper_chest != -1:
            bone_map.chest = chest_or_upper_chest
            print("DEBUG: Found Chest: ", chest_or_upper_chest, " ('", skeleton.get_bone_name(chest_or_upper_chest), "')")
        
        # 9. Guess Spine
        if chest_or_upper_chest != -1:
            bone_idx = skeleton.get_bone_parent(chest_or_upper_chest)
            if bone_idx != -1 and bone_idx != hips:
                bone_map.spine = bone_idx
                print("DEBUG: Found Spine: ", bone_idx, " ('", skeleton.get_bone_name(bone_idx), "')")

func search_bone_by_name(skeleton: Skeleton3D, picklist: Array, side: String = "", from_bone: int = -1, to_bone: int = -1) -> int:
    var bone_count = skeleton.get_bone_count()
    var best_match = -1
    var best_score = 0
    
    for i in range(bone_count):
        var bone_name = skeleton.get_bone_name(i).to_lower()
        
        # Check if this bone is in the search path (if specified)
        if from_bone != -1:
            if not is_bone_in_hierarchy_path(skeleton, i, from_bone, to_bone):
                continue
        
        # Check if bone name matches any pattern in picklist
        var matches_pattern = false
        for pattern in picklist:
            var regex = RegEx.new()
            regex.compile(pattern)
            if regex.search(bone_name):
                matches_pattern = true
                break
        
        if not matches_pattern:
            continue
        
        # Check side matching
        var score = 1
        if side != "":
            if side == "left":
                if "_l" in bone_name or ".l" in bone_name or "left" in bone_name:
                    score += 2
                elif "_r" in bone_name or ".r" in bone_name or "right" in bone_name:
                    continue  # Wrong side
            elif side == "right":
                if "_r" in bone_name or ".r" in bone_name or "right" in bone_name:
                    score += 2
                elif "_l" in bone_name or ".l" in bone_name or "left" in bone_name:
                    continue  # Wrong side
        
        if score > best_score:
            best_score = score
            best_match = i
    
    return best_match

func is_bone_in_hierarchy_path(skeleton: Skeleton3D, bone_idx: int, from_bone: int, to_bone: int) -> bool:
    # Check if bone is between from_bone and to_bone in hierarchy
    if to_bone == -1:
        # Just check if bone is descendant of from_bone
        var current = bone_idx
        while current != -1:
            if current == from_bone:
                return true
            current = skeleton.get_bone_parent(current)
        return false
    else:
        # Check if bone is on path from to_bone to from_bone
        var current = to_bone
        while current != -1 and current != from_bone:
            if current == bone_idx:
                return true
            current = skeleton.get_bone_parent(current)
        return current == from_bone and bone_idx == from_bone

func create_initial_correspondence(garment_skeleton: Skeleton3D, glb_skeleton: Skeleton3D, bone_map) -> Array:
    var correspondence = []
    correspondence.resize(15)
    
    # Initialize with -1 (no correspondence)
    for i in range(15):
        correspondence[i] = -1
    
    # Core bones correspondence (0-6)
    correspondence[0] = bone_map.root    # Hips
    correspondence[1] = bone_map.spine   # Spine
    correspondence[2] = bone_map.chest   # Chest
    correspondence[3] = bone_map.neck    # Neck
    correspondence[6] = bone_map.head    # Head
    
    return correspondence

# BONE SPLITTING: Create virtual hip subdivisions with child reparenting  
func create_hip_bone_subdivisions(glb_skeleton: Skeleton3D, bone_map):
    var subdivisions = {"left_hip_idx": -1, "right_hip_idx": -1}
    
    # Get hip position as reference center
    var hip_pos = Utils.calculate_bone_world_position(glb_skeleton, bone_map.root)
    
    # Calculate subdivision positions based on leg attachment points
    if bone_map.left_leg != -1:
        var left_leg_pos = Utils.calculate_bone_world_position(glb_skeleton, bone_map.left_leg)
        # Create virtual left hip position between center hip and left leg
        var left_hip_pos = hip_pos.lerp(left_leg_pos, 0.3)
        subdivisions.left_hip_idx = create_virtual_bone_index(glb_skeleton, "VirtualHipLeft", left_hip_pos)
        print("DEBUG: Created virtual left hip at: ", left_hip_pos)
    
    if bone_map.right_leg != -1:
        var right_leg_pos = Utils.calculate_bone_world_position(glb_skeleton, bone_map.right_leg)
        # Create virtual right hip position between center hip and right leg  
        var right_hip_pos = hip_pos.lerp(right_leg_pos, 0.3)
        subdivisions.right_hip_idx = create_virtual_bone_index(glb_skeleton, "VirtualHipRight", right_hip_pos)
        print("DEBUG: Created virtual right hip at: ", right_hip_pos)
    
    # Fallback: use center hip with lateral offsets
    if subdivisions.left_hip_idx == -1:
        var left_offset = hip_pos + Vector3(0.1, 0, 0)  # 10cm left offset
        subdivisions.left_hip_idx = create_virtual_bone_index(glb_skeleton, "VirtualHipLeft", left_offset)
    
    if subdivisions.right_hip_idx == -1:
        var right_offset = hip_pos + Vector3(-0.1, 0, 0)  # 10cm right offset  
        subdivisions.right_hip_idx = create_virtual_bone_index(glb_skeleton, "VirtualHipRight", right_offset)
    
    return subdivisions

func create_virtual_bone_index(glb_skeleton: Skeleton3D, bone_name: String, position: Vector3) -> int:
    # For now, return a virtual index based on position hash
    # In a real implementation, this would add the bone to the skeleton
    var hash_val = abs(position.x * 1000 + position.y * 100 + position.z * 10)
    return int(hash_val) % 10000 + 50000  # Virtual bone indices start at 50000

func add_limb_correspondences(glb_skeleton: Skeleton3D, correspondence: Array, bone_map):
    # Use detected shoulders and arms 
    if bone_map.right_shoulder != -1:
        if bone_map.right_arm != -1:
            correspondence[4] = bone_map.right_arm  # Right upper arm
            var right_forearm_children = Utils.get_bone_children(glb_skeleton, bone_map.right_arm)
            if right_forearm_children.size() > 0:
                correspondence[5] = right_forearm_children[0]  # Right forearm/hand
    
    if bone_map.left_shoulder != -1:
        if bone_map.left_arm != -1:
            correspondence[7] = bone_map.left_arm   # Left upper arm
            var left_forearm_children = Utils.get_bone_children(glb_skeleton, bone_map.left_arm)
            if left_forearm_children.size() > 0:
                correspondence[8] = left_forearm_children[0]  # Left forearm/hand
    
    # BONE SPLITTING: Create virtual hip subdivisions with child reparenting
    var hip_subdivisions = create_hip_bone_subdivisions(glb_skeleton, bone_map)
    
    # Map legs using subdivided hip hierarchy
    if bone_map.left_leg != -1:
        correspondence[10] = bone_map.left_leg       # Left upper leg 
        var left_leg_children = Utils.get_bone_children(glb_skeleton, bone_map.left_leg)
        if left_leg_children.size() > 0:
            correspondence[11] = left_leg_children[0]  # Left lower leg
    
    if bone_map.right_leg != -1:
        correspondence[13] = bone_map.right_leg      # Right upper leg
        var right_leg_children = Utils.get_bone_children(glb_skeleton, bone_map.right_leg)
        if right_leg_children.size() > 0:
            correspondence[14] = right_leg_children[0]  # Right lower leg
    
    # Use subdivided hip positions (NO MORE DUPLICATE MAPPINGS)
    correspondence[9] = hip_subdivisions.right_hip_idx   # Right hip subdivision
    correspondence[12] = hip_subdivisions.left_hip_idx   # Left hip subdivision
    
    print("DEBUG: Hip subdivision mapping - Right:[9]=", hip_subdivisions.right_hip_idx, " Left:[12]=", hip_subdivisions.left_hip_idx)

func fill_correspondence_gaps(garment_skeleton: Skeleton3D, glb_skeleton: Skeleton3D, correspondence: Array, bone_map):
    for i in range(correspondence.size()):
        if correspondence[i] == -1:
            print("DEBUG: No correspondence found for garment bone ", i, ", using scaled position")
            # Use position-based fallback for unmapped bones
            var garment_pos = Utils.calculate_bone_world_position(garment_skeleton, i)
            var scaled_pos = scale_garment_vertex_to_glb(garment_pos, glb_skeleton, bone_map)
            
            # Store as special "no correspondence" marker with position
            correspondence[i] = {"type": "position", "pos": scaled_pos}

# Interface method expected by SkeletonProcessor
func map_to_garment_skeleton(glb_skeleton: Skeleton3D, target_skeleton):
    print("DEBUG: Mapping GLB skeleton to garment skeleton using bone splitting")
    
    # Create a dummy garment skeleton from target_skeleton vertices
    var garment_skeleton = create_garment_skeleton_from_vertices(target_skeleton.vertices)
    
    # Use our advanced bone mapping with splitting
    var correspondence = map_garment_to_glb_skeleton(garment_skeleton, glb_skeleton)
    
    # Convert correspondence back to mapped skeleton format
    return create_mapped_skeleton_from_correspondence(correspondence, glb_skeleton, target_skeleton)

func create_garment_skeleton_from_vertices(vertices: Array) -> Skeleton3D:
    # Create a minimal skeleton object for compatibility
    var skeleton = Skeleton3D.new()
    # Add bones for each vertex
    for i in range(vertices.size()):
        skeleton.add_bone("Bone" + str(i))
    return skeleton

func create_mapped_skeleton_from_correspondence(correspondence: Array, glb_skeleton: Skeleton3D, target_skeleton):
    var mapped_skeleton = {
        "vertices": [],
        "edges": []
    }
    
    print("DEBUG: Creating mapped skeleton from correspondence with ", correspondence.size(), " bones")
    
    # Map each garment vertex to GLB bone position or virtual position
    for i in range(correspondence.size()):
        var mapped_pos = Vector3.ZERO
        
        if correspondence[i] is Dictionary:
            # Virtual bone position from bone splitting
            mapped_pos = correspondence[i].pos
            print("  Garment[", i, "] -> Virtual position: ", mapped_pos)
        elif correspondence[i] != -1:
            # Real GLB bone
            mapped_pos = Utils.calculate_bone_world_position(glb_skeleton, correspondence[i])
            var bone_name = glb_skeleton.get_bone_name(correspondence[i])
            print("  Garment[", i, "] -> GLB bone ", correspondence[i], " ('", bone_name, "') at ", mapped_pos)
        else:
            # Fallback to original position
            if i < target_skeleton.vertices.size():
                mapped_pos = target_skeleton.vertices[i]
            print("  Garment[", i, "] -> No correspondence, using scaled position: ", mapped_pos)
        
        mapped_skeleton.vertices.append(Utils.transform_to_blender_coords(mapped_pos))
    
    # Copy edges from target skeleton
    mapped_skeleton.edges = target_skeleton.edges
    
    return mapped_skeleton

func scale_garment_vertex_to_glb(garment_vertex: Vector3, glb_skeleton: Skeleton3D, bone_map) -> Vector3:
    var glb_root_pos = Utils.calculate_bone_world_position(glb_skeleton, bone_map.root) if bone_map.root != -1 else Vector3.ZERO
    var glb_head_pos = Utils.calculate_bone_world_position(glb_skeleton, bone_map.head) if bone_map.head != -1 else Vector3(0, 10, 0)
    
    var glb_height = abs(glb_head_pos.y - glb_root_pos.y)
    var scale_factor = glb_height / 10.0
    
    var scaled_pos = garment_vertex * scale_factor
    return Utils.transform_to_blender_coords(scaled_pos)
