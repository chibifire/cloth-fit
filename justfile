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
    cd {{build_dir}} && cmake -DCMAKE_BUILD_TYPE={{build_type}} -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..

compile: configure
    cd {{build_dir}} && make -j

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
