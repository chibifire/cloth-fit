set(SOURCES
	MatrixIO.cpp
	MatrixIO.hpp
	MshReader.cpp
	MshReader.hpp
	OBJReader.cpp
	OBJReader.hpp
	OBJWriter.cpp
	OBJWriter.hpp
	YamlToJson.cpp
	YamlToJson.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(polyfem PRIVATE ${SOURCES})
