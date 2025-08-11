class_name Utils
extends RefCounted

static func transform_to_blender_coords(vertex: Vector3) -> Vector3:
    # Transform from lying flat (Z = body length) to upright (Y = height)
    # Original: X=width, Y=thickness, Z=length (lying flat)
    # Target: X=width, Y=height, Z=depth (standing upright)
    return vertex

static func calculate_bone_world_position(skeleton: Skeleton3D, bone_idx: int) -> Vector3:
    # Use global pose instead of rest pose for accurate current position
    var global_transform = skeleton.get_bone_global_pose(bone_idx)
    return global_transform.origin

static func get_bone_children(skeleton: Skeleton3D, bone_idx: int) -> Array:
    var children = []
    var bone_count = skeleton.get_bone_count()
    
    for i in range(bone_count):
        if skeleton.get_bone_parent(i) == bone_idx:
            children.append(i)
    
    return children

static func calculate_distance_from_root(skeleton: Skeleton3D, bone_idx: int) -> int:
    var distance = 0
    var current_bone = bone_idx
    
    while skeleton.get_bone_parent(current_bone) != -1:
        distance += 1
        current_bone = skeleton.get_bone_parent(current_bone)
    
    return distance

static func find_skeletons(node: Node) -> Array:
    var skeletons = []
    if node is Skeleton3D:
        skeletons.append(node)
    for child in node.get_children():
        skeletons.append_array(find_skeletons(child))
    return skeletons

static func find_mesh_instances(node: Node) -> Array:
    var instances = []
    if node is MeshInstance3D or node.get_class() == "ImporterMeshInstance3D":
        instances.append(node)
    for child in node.get_children():
        instances.append_array(find_mesh_instances(child))
    return instances
