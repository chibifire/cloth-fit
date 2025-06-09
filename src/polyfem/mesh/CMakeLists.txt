set(SOURCES
	Mesh.cpp
	Mesh.hpp
	MeshUtils.cpp
	MeshUtils.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(polyfem PRIVATE ${SOURCES})

add_subdirectory(mesh2D)
add_subdirectory(mesh3D)
