// glTF 2.0 reader implementation using cgltf library
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <vector>

#include <Eigen/Core>

namespace polyfem::io
{
	class GLTFReader
	{
	public:
		GLTFReader() = delete;

		/// @brief Read a mesh from a glTF file
		///
		/// Fills in vertex positions, normals and texture coordinates. Mesh may
		/// have faces of any number of degree. Respects glTF coordinate system
		/// convention (right-handed with Y-up).
		///
		/// @param[in] gltf_file_name  path to .gltf or .glb file
		/// @param[out] V             double matrix of vertex positions
		/// @param[out] TC            double matrix of texture coordinates
		/// @param[out] N             double matrix of corner normals #N by 3
		/// @param[out] F             #F list of face indices into vertex positions
		/// @param[out] FTC           #F list of face indices into vertex texture
		///                           coordinates
		/// @param[out] FN            #F list of face indices into vertex normals
		/// @param[out] L             list of polyline indices into vertex positions
		///
		/// @returns true on success, false on errors
		static bool read(
			const std::string gltf_file_name,
			std::vector<std::vector<double>> &V,
			std::vector<std::vector<double>> &TC,
			std::vector<std::vector<double>> &N,
			std::vector<std::vector<int>> &F,
			std::vector<std::vector<int>> &FTC,
			std::vector<std::vector<int>> &FN,
			std::vector<std::vector<int>> &L);

		/// @brief Just read V, F, and L from glTF file
		static bool read(
			const std::string gltf_file_name,
			std::vector<std::vector<double>> &V,
			std::vector<std::vector<int>> &F,
			std::vector<std::vector<int>> &L);

		/// @brief Eigen Wrappers of read_gltf.
		/// @returns These will return true only if the data is perfectly
		///          "rectangular": All faces are the same degree, all have the same
		///          number of textures/normals etc.
		static bool read(
			const std::string str,
			Eigen::MatrixXd &V,
			Eigen::MatrixXi &E,
			Eigen::MatrixXi &F);

	};
} // namespace polyfem::io
