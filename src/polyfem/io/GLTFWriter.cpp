// glTF 2.0 writer implementation using tinygltf library
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "GLTFWriter.hpp"

#define TINYGLTF_NO_STB_IMAGE_WRITE
#define TINYGLTF_NO_STB_IMAGE
#define TINYGLTF_NO_EXTERNAL_IMAGE
#include <tiny_gltf.h>

#include <polyfem/utils/Logger.hpp>

#include <fstream>
#include <vector>
#include <cstring>

namespace polyfem::io
{
	bool GLTFWriter::write(const std::string &path, const Eigen::MatrixXd &v, const Eigen::MatrixXi &e, const Eigen::MatrixXi &f)
	{
		tinygltf::Model model;
		tinygltf::TinyGLTF writer;

		// Create a scene
		tinygltf::Scene scene;
		scene.nodes.push_back(0); // Reference to the main node
		model.scenes.push_back(scene);
		model.defaultScene = 0;

		// Create a node
		tinygltf::Node node;
		node.mesh = 0;
		model.nodes.push_back(node);

		// Create a mesh
		tinygltf::Mesh mesh;
		mesh.name = "mesh";

		// Create vertex buffer data
		std::vector<float> vertex_data;
		vertex_data.reserve(v.rows() * 3);
		
		for (int i = 0; i < v.rows(); ++i)
		{
			vertex_data.push_back(static_cast<float>(v(i, 0)));
			vertex_data.push_back(static_cast<float>(v(i, 1)));
			if (v.cols() >= 3)
				vertex_data.push_back(static_cast<float>(v(i, 2)));
			else
				vertex_data.push_back(0.0f);
		}

		// Create buffer for vertex data
		tinygltf::Buffer vertex_buffer;
		vertex_buffer.data.resize(vertex_data.size() * sizeof(float));
		std::memcpy(vertex_buffer.data.data(), vertex_data.data(), vertex_buffer.data.size());
		model.buffers.push_back(vertex_buffer);

		// Create buffer view for vertex positions
		tinygltf::BufferView vertex_buffer_view;
		vertex_buffer_view.buffer = 0;
		vertex_buffer_view.byteOffset = 0;
		vertex_buffer_view.byteLength = vertex_buffer.data.size();
		vertex_buffer_view.target = TINYGLTF_TARGET_ARRAY_BUFFER;
		model.bufferViews.push_back(vertex_buffer_view);

		// Create accessor for vertex positions
		tinygltf::Accessor vertex_accessor;
		vertex_accessor.bufferView = 0;
		vertex_accessor.byteOffset = 0;
		vertex_accessor.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
		vertex_accessor.count = v.rows();
		vertex_accessor.type = TINYGLTF_TYPE_VEC3;
		vertex_accessor.minValues = {static_cast<double>(v.col(0).minCoeff()), 
									 static_cast<double>(v.col(1).minCoeff()), 
									 v.cols() >= 3 ? static_cast<double>(v.col(2).minCoeff()) : 0.0};
		vertex_accessor.maxValues = {static_cast<double>(v.col(0).maxCoeff()), 
									 static_cast<double>(v.col(1).maxCoeff()), 
									 v.cols() >= 3 ? static_cast<double>(v.col(2).maxCoeff()) : 0.0};
		model.accessors.push_back(vertex_accessor);

		// Handle faces if provided
		if (f.rows() > 0)
		{
			// Create index buffer data
			std::vector<uint32_t> index_data;
			index_data.reserve(f.rows() * f.cols());
			
			for (int i = 0; i < f.rows(); ++i)
			{
				for (int j = 0; j < f.cols(); ++j)
				{
					index_data.push_back(static_cast<uint32_t>(f(i, j)));
				}
			}

			// Create buffer for index data
			tinygltf::Buffer index_buffer;
			index_buffer.data.resize(index_data.size() * sizeof(uint32_t));
			std::memcpy(index_buffer.data.data(), index_data.data(), index_buffer.data.size());
			model.buffers.push_back(index_buffer);

			// Create buffer view for indices
			tinygltf::BufferView index_buffer_view;
			index_buffer_view.buffer = 1;
			index_buffer_view.byteOffset = 0;
			index_buffer_view.byteLength = index_buffer.data.size();
			index_buffer_view.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;
			model.bufferViews.push_back(index_buffer_view);

			// Create accessor for indices
			tinygltf::Accessor index_accessor;
			index_accessor.bufferView = 1;
			index_accessor.byteOffset = 0;
			index_accessor.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
			index_accessor.count = index_data.size();
			index_accessor.type = TINYGLTF_TYPE_SCALAR;
			model.accessors.push_back(index_accessor);

			// Create primitive
			tinygltf::Primitive primitive;
			primitive.attributes["POSITION"] = 0; // Position accessor index
			primitive.indices = 1; // Index accessor index
			primitive.mode = TINYGLTF_MODE_TRIANGLES;
			mesh.primitives.push_back(primitive);
		}

		// Handle edges if provided (and no faces)
		if (e.rows() > 0 && f.rows() == 0)
		{
			// Create index buffer data for lines
			std::vector<uint32_t> line_index_data;
			line_index_data.reserve(e.rows() * 2);
			
			for (int i = 0; i < e.rows(); ++i)
			{
				line_index_data.push_back(static_cast<uint32_t>(e(i, 0)));
				line_index_data.push_back(static_cast<uint32_t>(e(i, 1)));
			}

			// Create buffer for line index data
			tinygltf::Buffer line_index_buffer;
			line_index_buffer.data.resize(line_index_data.size() * sizeof(uint32_t));
			std::memcpy(line_index_buffer.data.data(), line_index_data.data(), line_index_buffer.data.size());
			model.buffers.push_back(line_index_buffer);

			// Create buffer view for line indices
			tinygltf::BufferView line_index_buffer_view;
			line_index_buffer_view.buffer = model.buffers.size() - 1;
			line_index_buffer_view.byteOffset = 0;
			line_index_buffer_view.byteLength = line_index_buffer.data.size();
			line_index_buffer_view.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;
			model.bufferViews.push_back(line_index_buffer_view);

			// Create accessor for line indices
			tinygltf::Accessor line_index_accessor;
			line_index_accessor.bufferView = model.bufferViews.size() - 1;
			line_index_accessor.byteOffset = 0;
			line_index_accessor.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
			line_index_accessor.count = line_index_data.size();
			line_index_accessor.type = TINYGLTF_TYPE_SCALAR;
			model.accessors.push_back(line_index_accessor);

			// Create primitive for lines
			tinygltf::Primitive line_primitive;
			line_primitive.attributes["POSITION"] = 0; // Position accessor index
			line_primitive.indices = model.accessors.size() - 1; // Line index accessor index
                        line_primitive.mode = TINYGLTF_MODE_LINE;
			mesh.primitives.push_back(line_primitive);
		}

		// If no faces or edges, create a primitive with just points
		if (f.rows() == 0 && e.rows() == 0)
		{
			tinygltf::Primitive point_primitive;
			point_primitive.attributes["POSITION"] = 0; // Position accessor index
			point_primitive.mode = TINYGLTF_MODE_POINTS;
			mesh.primitives.push_back(point_primitive);
		}

		model.meshes.push_back(mesh);

		// Write the file
		bool ret = false;
		if (path.substr(path.find_last_of(".") + 1) == "glb")
		{
			ret = writer.WriteGltfSceneToFile(&model, path, false, false, true, false);
		}
		else
		{
			ret = writer.WriteGltfSceneToFile(&model, path, false, false, false, false);
		}

		if (!ret)
		{
			logger().error("GLTFWriter::write: Failed to write glTF file: {}", path);
			return false;
		}

		return true;
	}

} // namespace polyfem::io
