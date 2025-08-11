set(SOURCES
	MatrixIO.cpp
	MatrixIO.hpp
	OBJReader.cpp
	OBJReader.hpp
	OBJWriter.cpp
	OBJWriter.hpp
	GLTFReader.cpp
	GLTFReader.hpp
	GLTFWriter.cpp
	GLTFWriter.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(polyfem PRIVATE ${SOURCES})
