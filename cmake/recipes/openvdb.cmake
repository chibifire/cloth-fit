# Polyfem Solvers (https://github.com/Huangzizhou/openvdb)
# License: MIT

if(TARGET openvdb::openvdb)
    return()
endif()

message(STATUS "Third-party: creating target 'openvdb::openvdb'")

option(OPENVDB_CORE_SHARED "Build dynamically linked version of the core library." OFF)
option(OPENVDB_ENABLE_UNINSTALL "Adds a CMake uninstall target." OFF)

include(CPM)
CPMAddPackage("gh:Huangzizhou/openvdb#f9e84246f17b0870451dfd27f97bc8d1a92ec8f4")
