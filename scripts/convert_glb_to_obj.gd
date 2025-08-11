extends SceneTree

const GLTFLoader = preload("GLTFLoader.gd")
const MeshExporter = preload("MeshExporter.gd")
const SkeletonProcessor = preload("SkeletonProcessor.gd")
const WeightExtractor = preload("WeightExtractor.gd")

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
    
    var loader = GLTFLoader.new()
    var gltf_data = loader.load_glb(glb_file)
    
    if not gltf_data.success:
        print("Error loading GLB file: ", gltf_data.error)
        temp_node.queue_free()
        quit(1)
        return
    
    temp_node.add_child(gltf_data.scene)
    
    # Calculate avatar transformation to match garment skeleton
    var avatar_transform = {
        "scale_factor": 1.0,
        "hip_offset": Vector3.ZERO
    }
    if target_avatar != "":
        var skeleton_processor = SkeletonProcessor.new()
        avatar_transform = skeleton_processor.calculate_avatar_to_garment_transform(gltf_data.scene)
    
    var mesh_exporter = MeshExporter.new()
    mesh_exporter.export_meshes_as_obj(gltf_data.scene, output_name + ".obj", avatar_transform)
    mesh_exporter.create_mtl_file(output_name + ".mtl")
    
    var skeleton_processor = SkeletonProcessor.new()
    skeleton_processor.extract_armature_skeleton(gltf_data.scene, output_name + "_skeleton.obj", avatar_transform)
    
    if target_avatar != "":
        skeleton_processor.create_overlay_skeleton(gltf_data.scene, output_name + "_" + target_avatar.to_lower() + "_overlay.obj", target_avatar, garment_type)
    
    var weight_extractor = WeightExtractor.new()
    weight_extractor.extract_skinning_weights(gltf_data.scene, gltf_data.gltf_state, "skin.txt")
    
    print("Converted ", glb_file, " to:")
    print("  - ", output_name, ".obj (mesh)")
    print("  - ", output_name, ".mtl (materials)")
    print("  - ", output_name, "_skeleton.obj (skeleton)")
    if target_avatar != "":
        print("  - ", output_name, "_", target_avatar.to_lower(), "_overlay.obj (overlay skeleton)")
    print("  - skin.txt (skinning weights)")
    
    temp_node.queue_free()
    quit(0)
