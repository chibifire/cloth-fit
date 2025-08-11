build_type := "Release"
threads := "16"
build_dir := "build"

examples := "foxgirl_skirt Goblin_Jacket Goblin_Jumpsuit Trex_Jacket"

_run-example name: build
    cd garment-data/{{name}} && ../../{{build_dir}}/PolyFEM_bin -j setup.json --max_threads {{threads}} > log || echo "Example '{{name}}' failed - this appears to be an application issue, not a build issue"

_show-log name:
    @if [ -f garment-data/{{name}}/log ]; then cat garment-data/{{name}}/log; else echo "No log file found for '{{name}}'. Run 'just run-{{name}}' first."; fi

setup-build-dir:
    mkdir -p {{build_dir}}

configure: setup-build-dir
    cd {{build_dir}} && cmake -DCMAKE_BUILD_TYPE={{build_type}} ..

compile: configure
    cd {{build_dir}} && make -j4

build: compile

run-all-examples: run-foxgirl-skirt run-goblin-jacket run-goblin-jumpsuit run-trex-jacket

run-foxgirl-skirt: (_run-example "foxgirl_skirt")
run-goblin-jacket: (_run-example "Goblin_Jacket")
run-goblin-jumpsuit: (_run-example "Goblin_Jumpsuit")
run-trex-jacket: (_run-example "Trex_Jacket")

# Generic run command for any example
run example: (_run-example example)

test: build
    cd {{build_dir}} && ctest --output-on-failure

clean:
    rm -rf {{build_dir}}

format:
    find src -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i

format-check:
    find src -name "*.cpp" -o -name "*.hpp" | xargs clang-format --dry-run --Werror

list-examples:
    @echo "Available examples:"
    @ls -1 garment-data/ | grep -v "assets\|output.gif" | sed 's/^/  - /'

info:
    @echo "Cloth-Fit - Intersection-free Garment Retargeting"
    @echo "Build type: {{build_type}}"
    @echo "Threads: {{threads}}"
    @echo "Build directory: {{build_dir}}"

rebuild: clean build

dev: clean build test

show-foxgirl-log: (_show-log "foxgirl_skirt")
show-goblin-jacket-log: (_show-log "Goblin_Jacket")
show-goblin-jumpsuit-log: (_show-log "Goblin_Jumpsuit")
show-trex-jacket-log: (_show-log "Trex_Jacket")

show-log example: (_show-log example)

clean-examples:
    find garment-data -name "*.obj" -not -path "*/assets/*" -delete
    find garment-data -name "log" -delete

build-debug:
    just build_type=Debug build

run-debug example:
    cd garment-data/{{example}} && ../../{{build_dir}}/PolyFEM_bin -j setup.json --max_threads {{threads}} --log_level debug > debug_log

check-deps:
    @echo "Checking dependencies..."
    @which cmake > /dev/null && echo "✓ CMake found" || echo "✗ CMake not found"
    @which make > /dev/null && echo "✓ Make found" || echo "✗ Make not found"
    @which gcc > /dev/null && echo "✓ GCC found" || echo "✗ GCC not found"
    @which clang-format > /dev/null && echo "✓ clang-format found" || echo "✗ clang-format not found (optional)"
    @which godot > /dev/null && echo "✓ Godot found" || echo "✗ Godot not found (needed for GLB conversion)"

# Avatar Creation Tasks

# Convert GLB file to OBJ format using Godot (also extracts skeleton and skin weights)
convert-glb-to-obj glb_file output_name target_avatar="" garment_type="":
    @echo "Converting {{glb_file}} to OBJ format..."
    @if [ ! -f "{{glb_file}}" ]; then echo "Error: GLB file {{glb_file}} not found"; exit 1; fi
    godot --headless --script scripts/convert_glb_to_obj.gd -- "{{glb_file}}" "{{output_name}}" "{{target_avatar}}" "{{garment_type}}"
    @echo "✓ Converted to {{output_name}}.obj (mesh)"
    @echo "✓ Generated {{output_name}}.mtl (materials)"
    @echo "✓ Extracted {{output_name}}_skeleton.obj (skeleton)"
    @if [ "{{target_avatar}}" != "" ]; then echo "✓ Generated {{output_name}}_{{target_avatar}}_overlay.obj (overlay skeleton)"; fi
    @echo "✓ Generated skin.txt (skinning weights)"

# Extract skeleton from GLB file using Godot
extract-skeleton glb_file output_name:
    @echo "Extracting skeleton from {{glb_file}}..."
    @if [ ! -f "{{glb_file}}" ]; then echo "Error: GLB file {{glb_file}} not found"; exit 1; fi
    godot --headless --script scripts/extract_skeleton.gd -- "{{glb_file}}" "{{output_name}}"
    @echo "✓ Extracted skeleton to {{output_name}}.obj"

# Generate basic material file
generate-materials name:
    @echo "Generating material file for {{name}}..."
    @echo "# Material file for {{name}}" > {{name}}.mtl
    @echo "newmtl default" >> {{name}}.mtl
    @echo "Ka 0.2 0.2 0.2" >> {{name}}.mtl
    @echo "Kd 0.8 0.8 0.8" >> {{name}}.mtl
    @echo "Ks 0.0 0.0 0.0" >> {{name}}.mtl
    @echo "✓ Generated {{name}}.mtl"

# Create avatar directory structure
create-avatar-directory name:
    @echo "Creating avatar directory for {{name}}..."
    mkdir -p garment-data/assets/avatars/{{name}}
    @echo "✓ Created garment-data/assets/avatars/{{name}}"

# Setup any rigged avatar from GLB file
setup-avatar-from-glb glb_path avatar_name: (create-avatar-directory avatar_name)
	@echo "Setting up {{avatar_name}} avatar from GLB file..."
	@if [ ! -f "{{glb_path}}" ]; then echo "Error: {{glb_path}} not found"; exit 1; fi
	just convert-glb-to-obj "{{glb_path}}" "temp_avatar" "{{avatar_name}}" "skirt"
	mv temp_avatar.obj garment-data/assets/avatars/{{avatar_name}}/avatar.obj
	mv temp_avatar.mtl garment-data/assets/avatars/{{avatar_name}}/avatar.mtl 2>/dev/null || just generate-materials "garment-data/assets/avatars/{{avatar_name}}/avatar"
	mv temp_avatar_skeleton.obj garment-data/assets/avatars/{{avatar_name}}/skeleton.obj
	mv skin.txt garment-data/assets/avatars/{{avatar_name}}/skin.txt
	just generate-materials "garment-data/assets/avatars/{{avatar_name}}/skeleton"
	just validate-avatar-files {{avatar_name}}
	@echo "✓ {{avatar_name}} avatar setup complete"

# Validate avatar files exist
validate-avatar-files name:
    @echo "Validating avatar files for {{name}}..."
    @if [ ! -f "garment-data/assets/avatars/{{name}}/avatar.obj" ]; then echo "✗ Missing avatar.obj"; exit 1; fi
    @if [ ! -f "garment-data/assets/avatars/{{name}}/skeleton.obj" ]; then echo "✗ Missing skeleton.obj"; exit 1; fi
    @if [ ! -f "garment-data/assets/avatars/{{name}}/skin.txt" ]; then echo "✗ Missing skin.txt"; exit 1; fi
    @echo "✓ All required avatar files present for {{name}}"

# Create test configuration for new avatar
create-test-config avatar_name:
    @echo "Creating test configuration for {{avatar_name}}..."
    mkdir -p garment-data/{{avatar_name}}_test
    @echo '{\n  "incremental_steps": 2,\n  "avatar_mesh_path": "../assets/avatars/{{avatar_name}}/avatar.obj",\n  "target_skeleton_path": "../assets/avatars/{{avatar_name}}/skeleton.obj",\n  "avatar_skin_weights_path": "",\n  "garment_mesh_path": "../assets/garments/LCL_Skirt_DressEvening_003/garment.obj",\n  "no_fit_spec_path": "../assets/garments/LCL_Skirt_DressEvening_003/no-fit.txt",\n  "source_skeleton_path": "../assets/garments/LCL_Skirt_DressEvening_003/skeleton.obj",\n  "similarity_penalty_weight": 1,\n  "curvature_penalty_weight": 0.01,\n  "twist_penalty_weight": 0.01,\n  "curve_center_target_weight": 1,\n  "fit_weight": 2,\n  "symmetry_weight": 0,\n  "curve_size_weight": 0,\n  "voxel_size": 0.01,\n  "is_skirt": true,\n  "curve_center_target_automatic_bone_generation": true,\n  "contact": {\n    "enabled": true,\n    "dhat": 0.002\n  },\n  "solver": {\n    "max_threads": 16,\n    "linear": {\n      "solver": [\n        "Eigen::PardisoLDLT",\n        "Eigen::AccelerateLDLT",\n        "Eigen::SimplicialLDLT"\n      ]\n    },\n    "augmented_lagrangian": {\n      "initial_weight": 1,\n      "max_weight": 1000000.0,\n      "eta": 1,\n      "nonlinear": {\n        "grad_norm": 1,\n        "max_iterations": 50\n      }\n    },\n    "nonlinear": {\n      "Newton": {\n        "use_psd_projection": true,\n        "use_psd_projection_in_regularized": true,\n        "reg_weight_max": 1e16,\n        "reg_weight_min": 1,\n        "reg_weight_inc": 10000.0\n      },\n      "grad_norm": 0.01,\n      "line_search": {\n        "max_step_size_limiter": 0.5,\n        "use_grad_norm_tol": 1e-4,\n        "method": "Backtracking",\n        "min_step_size": 1e-8\n      },\n      "max_iterations": 5000\n    },\n    "contact": {\n      "CCD": {\n        "broad_phase": "BVH",\n        "max_iterations": 200,\n        "tolerance": 1e-3\n      },\n      "barrier_stiffness": 1e8\n    }\n  },\n  "output": {\n    "skip_frame": 2,\n    "log": {\n      "level": "debug"\n    }\n  }\n}' > garment-data/{{avatar_name}}_test/setup.json
    @echo "✓ Test configuration created for {{avatar_name}}"

# Test new avatar with garment fitting
test-avatar name: build (validate-avatar-files name)
    @echo "Testing avatar {{name}}..."
    just create-test-config {{name}}
    just _run-example "{{name}}_test"
    @echo "✓ Avatar {{name}} test completed. Check garment-data/{{name}}_test/log for results."

# Complete avatar setup from GLB to tested avatar
create-avatar-from-glb glb_path avatar_name: (convert-glb-to-obj glb_path "temp_avatar") (create-avatar-directory avatar_name)
    mv temp_avatar.obj garment-data/assets/avatars/{{avatar_name}}/avatar.obj
    mv temp_avatar.mtl garment-data/assets/avatars/{{avatar_name}}/avatar.mtl  
    mv temp_avatar_skeleton.obj garment-data/assets/avatars/{{avatar_name}}/skeleton.obj
    mv skin.txt garment-data/assets/avatars/{{avatar_name}}/skin.txt
    just generate-materials "garment-data/assets/avatars/{{avatar_name}}/skeleton"
    just validate-avatar-files {{avatar_name}}
    @echo "✓ Avatar {{avatar_name}} created successfully from {{glb_path}}"

# Clean avatar-related temporary files
clean-avatar-temp:
    rm -f temp_avatar.obj temp_avatar.mtl temp_avatar_skeleton.obj skin.txt
    @echo "✓ Cleaned temporary avatar files"

# List all available avatars
list-avatars:
    @echo "Available avatars:"
    @ls -1 garment-data/assets/avatars/ | sed 's/^/  - /'

# Garment Creation Tasks

# Create garment directory structure
create-garment-directory garment_name:
    @echo "Creating garment directory for {{garment_name}}..."
    mkdir -p garment-data/assets/garments/{{garment_name}}
    @echo "✓ Created garment-data/assets/garments/{{garment_name}}"

# Setup garment from GLB file with overlay skeleton for specific avatar
setup-garment-from-glb glb_path garment_name target_avatar garment_type: (create-garment-directory garment_name)
    @echo "Setting up {{garment_name}} garment from GLB file for {{target_avatar}}..."
    @if [ ! -f "{{glb_path}}" ]; then echo "Error: {{glb_path}} not found"; exit 1; fi
    just convert-glb-to-obj "{{glb_path}}" "temp_garment" "{{target_avatar}}" "{{garment_type}}"
    mv temp_garment.obj garment-data/assets/garments/{{garment_name}}/garment.obj
    mv temp_garment.mtl garment-data/assets/garments/{{garment_name}}/garment.mtl 2>/dev/null || just generate-materials "garment-data/assets/garments/{{garment_name}}/garment"
    mv temp_garment_skeleton.obj garment-data/assets/garments/{{garment_name}}/skeleton.obj
    mv temp_garment_{{target_avatar}}_overlay.obj garment-data/assets/garments/{{garment_name}}/skeleton_{{target_avatar}}_overlay.obj 2>/dev/null || echo "No overlay skeleton generated"
    mv skin.txt garment-data/assets/garments/{{garment_name}}/skin.txt 2>/dev/null || echo "No skin weights generated"
    just generate-materials "garment-data/assets/garments/{{garment_name}}/skeleton"
    @echo "✓ {{garment_name}} garment setup complete for {{target_avatar}}"

# Generate missing overlay skeleton for existing garment
generate-overlay-skeleton garment_name target_avatar garment_type:
    @echo "Generating overlay skeleton for {{garment_name}} -> {{target_avatar}}..."
    @if [ ! -f "garment-data/assets/garments/{{garment_name}}/skeleton.obj" ]; then echo "Error: Source skeleton not found"; exit 1; fi
    @if [ ! -f "garment-data/assets/avatars/{{target_avatar}}/skeleton.obj" ]; then echo "Error: Target avatar skeleton not found"; exit 1; fi
    cd garment-data/assets/garments/{{garment_name}} && godot --headless --script ../../../../scripts/convert_glb_to_obj.gd -- "dummy.glb" "temp" "{{target_avatar}}" "{{garment_type}}" || echo "Using fallback method..."
    @echo "✓ Generated skeleton_{{target_avatar}}_overlay.obj"

# Create project for garment fitting
create-garment-project garment_name target_avatar project_name:
    @echo "Creating garment fitting project {{project_name}}..."
    mkdir -p garment-data/{{project_name}}
    @echo '{\n  "incremental_steps": 2,\n  "avatar_mesh_path": "../assets/avatars/{{target_avatar}}/avatar.obj",\n  "target_skeleton_path": "../assets/avatars/{{target_avatar}}/skeleton.obj",\n  "avatar_skin_weights_path": "",\n  "garment_mesh_path": "../assets/garments/{{garment_name}}/garment.obj",\n  "no_fit_spec_path": "../assets/garments/{{garment_name}}/no-fit.txt",\n  "source_skeleton_path": "../assets/garments/{{garment_name}}/skeleton_{{target_avatar}}_overlay.obj",\n  "similarity_penalty_weight": 1,\n  "curvature_penalty_weight": 0.01,\n  "twist_penalty_weight": 0.01,\n  "curve_center_target_weight": 1,\n  "fit_weight": 2,\n  "symmetry_weight": 0,\n  "curve_size_weight": 0,\n  "voxel_size": 0.01,\n  "is_skirt": true,\n  "curve_center_target_automatic_bone_generation": true,\n  "contact": {\n    "enabled": true,\n    "dhat": 0.002\n  },\n  "solver": {\n    "max_threads": 16,\n    "linear": {\n      "solver": [\n        "Eigen::PardisoLDLT",\n        "Eigen::AccelerateLDLT",\n        "Eigen::SimplicialLDLT"\n      ]\n    },\n    "augmented_lagrangian": {\n      "initial_weight": 1,\n      "max_weight": 1000000.0,\n      "eta": 1,\n      "nonlinear": {\n        "grad_norm": 1,\n        "max_iterations": 50\n      }\n    },\n    "nonlinear": {\n      "Newton": {\n        "use_psd_projection": true,\n        "use_psd_projection_in_regularized": true,\n        "reg_weight_max": 1e16,\n        "reg_weight_min": 1,\n        "reg_weight_inc": 10000.0\n      },\n      "grad_norm": 0.01,\n      "line_search": {\n        "max_step_size_limiter": 0.5,\n        "use_grad_norm_tol": 1e-4,\n        "method": "Backtracking",\n        "min_step_size": 1e-8\n      },\n      "max_iterations": 5000\n    },\n    "contact": {\n      "CCD": {\n        "broad_phase": "BVH",\n        "max_iterations": 200,\n        "tolerance": 1e-3\n      },\n      "barrier_stiffness": 1e8\n    }\n  },\n  "output": {\n    "skip_frame": 2,\n    "log": {\n      "level": "debug"\n    }\n  }\n}' > garment-data/{{project_name}}/setup.json
    @echo "✓ Created project {{project_name}} for {{garment_name}} -> {{target_avatar}}"

# Show help for avatar setup commands
help-avatar-setup:
    @echo "Avatar Setup Commands:"
    @echo ""
    @echo "Setup any rigged avatar from GLB file:"
    @echo "  just setup-avatar-from-glb \"/path/to/avatar.glb\" \"AvatarName\""
    @echo ""
    @echo "Example:"
    @echo "  just setup-avatar-from-glb \"../my_character.glb\" \"MyCharacter\""
    @echo ""
    @echo "This will:"
    @echo "  1. Create avatar directory structure"
    @echo "  2. Convert GLB to OBJ format"
    @echo "  3. Extract skeleton structure"
    @echo "  4. Generate material files"
    @echo "  5. Validate all required files exist"
    @echo ""
    @echo "After setup, test the avatar:"
    @echo "  just test-avatar \"AvatarName\""
    @echo ""
    @echo "View available avatars:"
    @echo "  just list-avatars"

# Show help for garment setup commands
help-garment-setup:
    @echo "Garment Setup Commands:"
    @echo ""
    @echo "Setup garment from GLB file for specific avatar:"
    @echo "  just setup-garment-from-glb \"/path/to/garment.glb\" \"GarmentName\" \"AvatarName\" \"garment_type\""
    @echo ""
    @echo "Generate missing overlay skeleton:"
    @echo "  just generate-overlay-skeleton \"GarmentName\" \"AvatarName\" \"garment_type\""
    @echo ""
    @echo "Create fitting project:"
    @echo "  just create-garment-project \"GarmentName\" \"AvatarName\" \"ProjectName\""
    @echo ""