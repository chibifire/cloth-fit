#
# Copyright 2020 Adobe. All rights reserved.
# This file is licensed to you under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may obtain a copy
# of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under
# the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR REPRESENTATIONS
# OF ANY KIND, either express or implied. See the License for the specific language
# governing permissions and limitations under the License.
#

# tinygltf (https://github.com/syoyo/tinygltf)
# License: MIT

if(TARGET tinygltf::tinygltf)
    return()
endif()

message(STATUS "Third-party: creating target 'tinygltf::tinygltf'")

include(CPM)
CPMAddPackage(
    NAME tinygltf
    GITHUB_REPOSITORY syoyo/tinygltf
    GIT_TAG v2.8.21
    DOWNLOAD_ONLY TRUE
)

add_library(tinygltf_tinygltf INTERFACE)
add_library(tinygltf::tinygltf ALIAS tinygltf_tinygltf)

include(GNUInstallDirs)
target_include_directories(tinygltf_tinygltf SYSTEM INTERFACE
    $<BUILD_INTERFACE:${tinygltf_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# tinygltf is header-only but requires C++11
target_compile_features(tinygltf_tinygltf INTERFACE cxx_std_11)

# Define required preprocessor macros for tinygltf
target_compile_definitions(tinygltf_tinygltf INTERFACE
    TINYGLTF_NO_STB_IMAGE_WRITE
    TINYGLTF_NO_STB_IMAGE
    TINYGLTF_NO_EXTERNAL_IMAGE
)

# Install rules
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME tinygltf)
set_target_properties(tinygltf_tinygltf PROPERTIES EXPORT_NAME tinygltf)
install(FILES ${tinygltf_SOURCE_DIR}/tiny_gltf.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS tinygltf_tinygltf EXPORT tinygltf_Targets)
install(EXPORT tinygltf_Targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tinygltf NAMESPACE tinygltf::)
