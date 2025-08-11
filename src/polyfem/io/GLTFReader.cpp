// glTF 2.0 reader implementation using tinygltf library
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "GLTFReader.hpp"

#define TINYGLTF_IMPLEMENTATION
#define TINYGLTF_NO_STB_IMAGE_WRITE
#define TINYGLTF_NO_STB_IMAGE
#define TINYGLTF_NO_EXTERNAL_IMAGE
#include <tiny_gltf.h>

#include <cstdio>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <limits>

#include <igl/edges.h>
#include <igl/list_to_matrix.h>

#include <polyfem/utils/Logger.hpp>

namespace polyfem::io
{
	bool GLTFReader::read(
		const std::string gltf_file_name,
		std::vector<std::vector<double>> &V,
		std::vector<std::vector<double>> &TC,
		std::vector<std::vector<double>> &N,
		std::vector<std::vector<int>> &F,
		std::vector<std::vector<int>> &FTC,
		std::vector<std::vector<int>> &FN,
		std::vector<std::vector<int>> &L)
	{
		// Clear outputs
		V.clear();
		TC.clear();
		N.clear();
		F.clear();
		FTC.clear();
		FN.clear();
		L.clear();

		// Load glTF file
		tinygltf::Model model;
		tinygltf::TinyGLTF loader;
		std::string err;
		std::string warn;

		bool ret = false;
		if (gltf_file_name.substr(gltf_file_name.find_last_of(".") + 1) == "glb")
		{
			ret = loader.LoadBinaryFromFile(&model, &err, &warn, gltf_file_name);
		}
		else
		{
			ret = loader.LoadASCIIFromFile(&model, &err, &warn, gltf_file_name);
		}

		if (!warn.empty())
		{
			logger().warn("GLTFReader::read: {}", warn);
		}

		if (!err.empty())
		{
			logger().error("GLTFReader::read: {}", err);
		}

		if (!ret)
		{
			logger().error("GLTFReader::read: Failed to parse glTF file: {}", gltf_file_name);
			return false;
		}

		// Process each mesh in the glTF file
		for (size_t i = 0; i < model.meshes.size(); ++i)
		{
			const tinygltf::Mesh& mesh = model.meshes[i];
			
			// Process each primitive in the mesh
			for (size_t j = 0; j < mesh.primitives.size(); ++j)
			{
				const tinygltf::Primitive& primitive = mesh.primitives[j];
				
				// Only handle triangular primitives for now
				if (primitive.mode != TINYGLTF_MODE_TRIANGLES)
				{
					logger().warn("GLTFReader::read: skipping non-triangular primitive (mode: {})", primitive.mode);
					continue;
				}

				// Get vertex positions
				auto pos_it = primitive.attributes.find("POSITION");
				if (pos_it == primitive.attributes.end())
				{
					logger().warn("GLTFReader::read: primitive missing POSITION attribute");
					continue;
				}

				const tinygltf::Accessor& pos_accessor = model.accessors[pos_it->second];
				const tinygltf::BufferView& pos_bufferView = model.bufferViews[pos_accessor.bufferView];
				const tinygltf::Buffer& pos_buffer = model.buffers[pos_bufferView.buffer];

				// Extract vertex data
				size_t vertex_count = pos_accessor.count;
				size_t vertex_offset = V.size();

				if (pos_accessor.type != TINYGLTF_TYPE_VEC3 || pos_accessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT)
				{
					logger().error("GLTFReader::read: unsupported position format");
					continue;
				}

				const float* positions = reinterpret_cast<const float*>(
					&pos_buffer.data[pos_bufferView.byteOffset + pos_accessor.byteOffset]);

				for (size_t v = 0; v < vertex_count; ++v)
				{
					std::vector<double> vertex = {
						static_cast<double>(positions[v * 3 + 0]),
						static_cast<double>(positions[v * 3 + 1]),
						static_cast<double>(positions[v * 3 + 2])
					};
					V.push_back(vertex);
				}

				// Read normals if available
				auto normal_it = primitive.attributes.find("NORMAL");
				if (normal_it != primitive.attributes.end())
				{
					const tinygltf::Accessor& normal_accessor = model.accessors[normal_it->second];
					const tinygltf::BufferView& normal_bufferView = model.bufferViews[normal_accessor.bufferView];
					const tinygltf::Buffer& normal_buffer = model.buffers[normal_bufferView.buffer];

					if (normal_accessor.type == TINYGLTF_TYPE_VEC3 && 
						normal_accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT)
					{
						const float* normals = reinterpret_cast<const float*>(
							&normal_buffer.data[normal_bufferView.byteOffset + normal_accessor.byteOffset]);

						for (size_t v = 0; v < vertex_count; ++v)
						{
							std::vector<double> normal_vec = {
								static_cast<double>(normals[v * 3 + 0]),
								static_cast<double>(normals[v * 3 + 1]),
								static_cast<double>(normals[v * 3 + 2])
							};
							N.push_back(normal_vec);
						}
					}
				}

				// Read texture coordinates if available
				auto texcoord_it = primitive.attributes.find("TEXCOORD_0");
				if (texcoord_it != primitive.attributes.end())
				{
					const tinygltf::Accessor& texcoord_accessor = model.accessors[texcoord_it->second];
					const tinygltf::BufferView& texcoord_bufferView = model.bufferViews[texcoord_accessor.bufferView];
					const tinygltf::Buffer& texcoord_buffer = model.buffers[texcoord_bufferView.buffer];

					if (texcoord_accessor.type == TINYGLTF_TYPE_VEC2 && 
						texcoord_accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT)
					{
						const float* texcoords = reinterpret_cast<const float*>(
							&texcoord_buffer.data[texcoord_bufferView.byteOffset + texcoord_accessor.byteOffset]);

						for (size_t v = 0; v < vertex_count; ++v)
						{
							std::vector<double> tc_vec = {
								static_cast<double>(texcoords[v * 2 + 0]),
								static_cast<double>(texcoords[v * 2 + 1])
							};
							TC.push_back(tc_vec);
						}
					}
				}

				// Read indices
				if (primitive.indices >= 0)
				{
					const tinygltf::Accessor& index_accessor = model.accessors[primitive.indices];
					const tinygltf::BufferView& index_bufferView = model.bufferViews[index_accessor.bufferView];
					const tinygltf::Buffer& index_buffer = model.buffers[index_bufferView.buffer];

					size_t face_count = index_accessor.count / 3; // Assuming triangles
					
					if (index_accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT)
					{
						const uint16_t* indices = reinterpret_cast<const uint16_t*>(
							&index_buffer.data[index_bufferView.byteOffset + index_accessor.byteOffset]);

						for (size_t face = 0; face < face_count; ++face)
						{
							std::vector<int> face_indices = {
								static_cast<int>(indices[face * 3 + 0] + vertex_offset),
								static_cast<int>(indices[face * 3 + 1] + vertex_offset),
								static_cast<int>(indices[face * 3 + 2] + vertex_offset)
							};
							F.push_back(face_indices);
						}
					}
					else if (index_accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT)
					{
						const uint32_t* indices = reinterpret_cast<const uint32_t*>(
							&index_buffer.data[index_bufferView.byteOffset + index_accessor.byteOffset]);

						for (size_t face = 0; face < face_count; ++face)
						{
							std::vector<int> face_indices = {
								static_cast<int>(indices[face * 3 + 0] + vertex_offset),
								static_cast<int>(indices[face * 3 + 1] + vertex_offset),
								static_cast<int>(indices[face * 3 + 2] + vertex_offset)
							};
							F.push_back(face_indices);
						}
					}
					else if (index_accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE)
					{
						const uint8_t* indices = reinterpret_cast<const uint8_t*>(
							&index_buffer.data[index_bufferView.byteOffset + index_accessor.byteOffset]);

						for (size_t face = 0; face < face_count; ++face)
						{
							std::vector<int> face_indices = {
								static_cast<int>(indices[face * 3 + 0] + vertex_offset),
								static_cast<int>(indices[face * 3 + 1] + vertex_offset),
								static_cast<int>(indices[face * 3 + 2] + vertex_offset)
							};
							F.push_back(face_indices);
						}
					}

					// For texture coordinates and normals, use same indices for now
					// This assumes vertex attributes are aligned (typical for most glTF files)
					for (size_t face = 0; face < face_count; ++face)
					{
						if (!TC.empty())
						{
							FTC.push_back(F[F.size() - face_count + face]);
						}
						else
						{
							FTC.push_back(std::vector<int>());
						}
						
						if (!N.empty())
						{
							FN.push_back(F[F.size() - face_count + face]);
						}
						else
						{
							FN.push_back(std::vector<int>());
						}
					}
				}
			}
		}

		return true;
	}

	bool GLTFReader::read(
		const std::string gltf_file_name,
		std::vector<std::vector<double>> &V,
		std::vector<std::vector<int>> &F,
		std::vector<std::vector<int>> &L)
	{
		std::vector<std::vector<double>> TC, N;
		std::vector<std::vector<int>> FTC, FN;
		return read(gltf_file_name, V, TC, N, F, FTC, FN, L);
	}

	bool GLTFReader::read(
		const std::string str,
		Eigen::MatrixXd &V,
		Eigen::MatrixXi &E,
		Eigen::MatrixXi &F)
	{
		std::vector<std::vector<double>> vV, vTC, vN;
		std::vector<std::vector<int>> vF, vFTC, vFN, vL;
		bool success = read(str, vV, vTC, vN, vF, vFTC, vFN, vL);
		if (!success)
		{
			return false;
		}

		bool V_rect = igl::list_to_matrix(vV, V);
		if (!V_rect)
		{
			logger().error("GLTFReader::read: vertices not rectangular matrix!");
			return false;
		}

		bool F_rect = igl::list_to_matrix(vF, F);
		if (!F_rect)
		{
			logger().error("GLTFReader::read: faces not rectangular matrix!");
			return false;
		}

		// Convert polylines to edges
		std::vector<std::vector<int>> vE;
		for (const std::vector<int> &polyline : vL)
		{
			for (size_t i = 1; i < polyline.size(); i++)
			{
				vE.push_back({{polyline[i - 1], polyline[i]}});
			}
		}

		bool E_rect = igl::list_to_matrix(vE, E);
		if (!E_rect)
		{
			logger().error("GLTFReader::read: edges not rectangular matrix!");
			return false;
		}

		return true;
	}

} // namespace polyfem::io
