// glTF 2.0 writer implementation using tinygltf library
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

#include <Eigen/Core>

namespace polyfem::io
{
	class GLTFWriter
	{
	public:
		GLTFWriter() = delete;

		/// @brief Write mesh data to glTF format
		/// @param[in] path Output file path (.gltf or .glb)
		/// @param[in] v Vertex positions matrix
		/// @param[in] e Edge indices matrix
		/// @param[in] f Face indices matrix
		/// @returns true on success, false on error
		static bool write(
			const std::string &path,
			const Eigen::MatrixXd &v,
			const Eigen::MatrixXi &e,
			const Eigen::MatrixXi &f);

		/// @brief Write mesh data with automatic type detection
		/// @param[in] path Output file path
		/// @param[in] v Vertex positions matrix
		/// @param[in] e_or_f Either edge indices (if 2 columns) or face indices (if 3+ columns)
		/// @returns true on success, false on error
		static bool write(
			const std::string &path,
			const Eigen::MatrixXd &v,
			const Eigen::MatrixXi &e_or_f)
		{
			if (e_or_f.cols() == 2)
				return write(path, v, e_or_f, Eigen::MatrixXi());
			else
				return write(path, v, Eigen::MatrixXi(), e_or_f);
		}
	};
} // namespace polyfem::io
