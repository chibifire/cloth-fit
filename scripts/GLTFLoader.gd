class_name GLTFLoader
extends RefCounted

func load_glb(file_path: String) -> Dictionary:
    var gltf_document = GLTFDocument.new()
    var gltf_state = GLTFState.new()
    var error = gltf_document.append_from_file(file_path, gltf_state)
    
    if error != OK:
        return {
            "success": false,
            "error": error,
            "scene": null,
            "gltf_state": null
        }
    
    var scene = gltf_document.generate_scene(gltf_state)
    
    if not scene:
        return {
            "success": false,
            "error": "Could not generate scene from GLB file",
            "scene": null,
            "gltf_state": null
        }
    
    return {
        "success": true,
        "error": null,
        "scene": scene,
        "gltf_state": gltf_state
    }
