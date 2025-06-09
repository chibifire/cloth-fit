set(AUTOGEN
	auto_derivatives1.cpp
	auto_derivatives2.cpp
	auto_derivatives3.cpp
	auto_derivatives.hpp
)

set(SOURCES
	auto_tetrahedron.ipp
	auto_triangle.ipp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${AUTOGEN})

target_sources(polyfem PRIVATE ${SOURCES} ${AUTOGEN})
