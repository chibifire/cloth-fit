#include "optimize.hpp"

#include <polyfem/solver/GarmentNLProblem.hpp>
#include <polyfem/solver/forms/ContactForm.hpp>
#include <polyfem/io/OBJReader.hpp>
#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/io/MatrixIO.hpp>
#include <polyfem/mesh/MeshUtils.hpp>
#include <polyfem/utils/JSONUtils.hpp>
#include <polyfem/utils/par_for.hpp>
#include <polyfem/utils/StringUtils.hpp>
#include <polyfem/utils/Logger.hpp>
#include <polyfem/utils/MatrixUtils.hpp>

#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readOBJ.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/write_triangle_mesh.h>
#include <igl/writeOBJ.h>

#include <polysolve/linear/Solver.hpp>

#ifdef POLYFEM_WITH_PARAVIEWO
#include <paraviewo/ParaviewWriter.hpp>
#include <paraviewo/VTUWriter.hpp>
#endif

#include <ipc/ipc.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/utils/logger.hpp>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>

#include <jse/jse.h>

#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <regex>
#include <iomanip>
#include <sstream>

using namespace polyfem::solver;
using namespace polyfem::mesh;

namespace spdlog::level
{
	NLOHMANN_JSON_SERIALIZE_ENUM(
		spdlog::level::level_enum,
		{{spdlog::level::level_enum::trace, "trace"},
		 {spdlog::level::level_enum::debug, "debug"},
		 {spdlog::level::level_enum::info, "info"},
		 {spdlog::level::level_enum::warn, "warning"},
		 {spdlog::level::level_enum::err, "error"},
		 {spdlog::level::level_enum::critical, "critical"},
		 {spdlog::level::level_enum::off, "off"},
		 {spdlog::level::level_enum::trace, 0},
		 {spdlog::level::level_enum::debug, 1},
		 {spdlog::level::level_enum::info, 2},
		 {spdlog::level::level_enum::warn, 3},
		 {spdlog::level::level_enum::err, 3},
		 {spdlog::level::level_enum::critical, 4},
		 {spdlog::level::level_enum::off, 5}})
}
namespace polyfem {
    namespace {
        /// @brief Parse MTL file and return materials map
        std::map<std::string, MTLMaterial> parse_mtl_file(const std::string &mtl_path)
        {
            std::map<std::string, MTLMaterial> materials;

            if (!std::filesystem::exists(mtl_path))
            {
                logger().warn("MTL file does not exist: {}", mtl_path);
                return materials;
            }

            std::ifstream mtl_file(mtl_path);
            if (!mtl_file.is_open())
            {
                logger().error("Failed to open MTL file: {}", mtl_path);
                return materials;
            }

            MTLMaterial current_material;
            bool has_current_material = false;
            std::string line;

            while (std::getline(mtl_file, line))
            {
                // Remove leading/trailing whitespace
                line.erase(0, line.find_first_not_of(" \t"));
                line.erase(line.find_last_not_of(" \t") + 1);

                // Skip empty lines and comments
                if (line.empty() || line[0] == '#')
                    continue;

                std::istringstream iss(line);
                std::string token;
                iss >> token;

                if (token == "newmtl")
                {
                    // Save previous material if exists
                    if (has_current_material && !current_material.name.empty())
                    {
                        materials[current_material.name] = current_material;
                    }

                    // Start new material
                    current_material = MTLMaterial();
                    iss >> current_material.name;
                    has_current_material = true;
                }
                else if (has_current_material)
                {
                    if (token == "Ka")
                    {
                        iss >> current_material.Ka[0] >> current_material.Ka[1] >> current_material.Ka[2];
                    }
                    else if (token == "Kd")
                    {
                        iss >> current_material.Kd[0] >> current_material.Kd[1] >> current_material.Kd[2];
                    }
                    else if (token == "Ks")
                    {
                        iss >> current_material.Ks[0] >> current_material.Ks[1] >> current_material.Ks[2];
                    }
                    else if (token == "Ke")
                    {
                        iss >> current_material.Ke[0] >> current_material.Ke[1] >> current_material.Ke[2];
                    }
                    else if (token == "Ns")
                    {
                        iss >> current_material.Ns;
                    }
                    else if (token == "Ni")
                    {
                        iss >> current_material.Ni;
                    }
                    else if (token == "d")
                    {
                        iss >> current_material.d;
                    }
                    else if (token == "illum")
                    {
                        iss >> current_material.illum;
                    }
                    else if (token == "map_Kd")
                    {
                        // Read the rest of the line to handle filenames with spaces
                        std::getline(iss, current_material.map_Kd);
                        // Remove leading whitespace
                        current_material.map_Kd.erase(0, current_material.map_Kd.find_first_not_of(" \t"));
                    }
                    else if (token == "map_d")
                    {
                        std::getline(iss, current_material.map_d);
                        current_material.map_d.erase(0, current_material.map_d.find_first_not_of(" \t"));
                    }
                    else if (token == "map_Ks")
                    {
                        std::getline(iss, current_material.map_Ks);
                        current_material.map_Ks.erase(0, current_material.map_Ks.find_first_not_of(" \t"));
                    }
                    else if (token == "map_Ka")
                    {
                        std::getline(iss, current_material.map_Ka);
                        current_material.map_Ka.erase(0, current_material.map_Ka.find_first_not_of(" \t"));
                    }
                    else if (token == "bump")
                    {
                        std::getline(iss, current_material.bump);
                        current_material.bump.erase(0, current_material.bump.find_first_not_of(" \t"));
                    }
                }
            }

            // Save last material
            if (has_current_material && !current_material.name.empty())
            {
                materials[current_material.name] = current_material;
            }

            mtl_file.close();
            logger().debug("Parsed {} materials from MTL file: {}", materials.size(), mtl_path);
            return materials;
        }

        /// @brief Write MTL file with materials and copy referenced textures
        bool write_mtl_file(const std::string &dest_mtl_path,
                           const std::map<std::string, MTLMaterial> &materials,
                           const std::string &source_mtl_dir = "")
        {
            try
            {
                // Create destination directory if it doesn't exist
                std::filesystem::create_directories(std::filesystem::path(dest_mtl_path).parent_path());

                std::ofstream mtl_file(dest_mtl_path);
                if (!mtl_file.is_open())
                {
                    logger().error("Failed to create MTL file: {}", dest_mtl_path);
                    return false;
                }

                // Write header comment
                mtl_file << "# Generated MTL file by PolyFEM cloth-fit\n";
                mtl_file << "# www.polyfem.github.io\n\n";

                std::filesystem::path dest_dir = std::filesystem::path(dest_mtl_path).parent_path();

                for (const auto &[name, material] : materials)
                {
                    mtl_file << "newmtl " << material.name << "\n";
                    mtl_file << "Ns " << std::fixed << std::setprecision(6) << material.Ns << "\n";
                    mtl_file << "Ka " << material.Ka[0] << " " << material.Ka[1] << " " << material.Ka[2] << "\n";
                    mtl_file << "Kd " << material.Kd[0] << " " << material.Kd[1] << " " << material.Kd[2] << "\n";
                    mtl_file << "Ks " << material.Ks[0] << " " << material.Ks[1] << " " << material.Ks[2] << "\n";
                    mtl_file << "Ke " << material.Ke[0] << " " << material.Ke[1] << " " << material.Ke[2] << "\n";
                    mtl_file << "Ni " << material.Ni << "\n";
                    mtl_file << "d " << material.d << "\n";
                    mtl_file << "illum " << material.illum << "\n";

                    // Handle texture maps and copy texture files
                    auto copy_texture = [&](const std::string &texture_name, const std::string &mtl_directive) {
                        if (!texture_name.empty())
                        {
                            mtl_file << mtl_directive << " " << texture_name << "\n";

                            // Copy texture file if source directory is provided
                            if (!source_mtl_dir.empty())
                            {
                                std::filesystem::path source_texture_path = std::filesystem::path(source_mtl_dir) / texture_name;
                                std::filesystem::path dest_texture_path = dest_dir / texture_name;

                                if (std::filesystem::exists(source_texture_path))
                                {
                                    try
                                    {
                                        std::filesystem::copy_file(source_texture_path, dest_texture_path,
                                                                  std::filesystem::copy_options::overwrite_existing);
                                        logger().debug("Copied texture file: {} -> {}", source_texture_path.string(), dest_texture_path.string());
                                    }
                                    catch (const std::exception& e)
                                    {
                                        logger().warn("Failed to copy texture file {}: {}", source_texture_path.string(), e.what());
                                    }
                                }
                                else
                                {
                                    logger().warn("Referenced texture file not found: {}", source_texture_path.string());
                                }
                            }
                        }
                    };

                    copy_texture(material.map_Kd, "map_Kd");
                    copy_texture(material.map_d, "map_d");
                    copy_texture(material.map_Ks, "map_Ks");
                    copy_texture(material.map_Ka, "map_Ka");
                    copy_texture(material.bump, "bump");

                    mtl_file << "\n";
                }

                mtl_file.close();
                logger().debug("Wrote MTL file with {} materials: {}", materials.size(), dest_mtl_path);
                return true;
            }
            catch (const std::exception& e)
            {
                logger().error("Failed to write MTL file {}: {}", dest_mtl_path, e.what());
                return false;
            }
        }

        Eigen::Vector2d project_to_line(
            const Eigen::Vector3d &a,
            const Eigen::Vector3d &b,
            const Eigen::Vector3d &p)
        {
            Eigen::Vector3d s = b - a;
            double t = (p - a).dot(s) / s.squaredNorm();
            double d = (p - (a + t * s)).squaredNorm();

            return Eigen::Vector2d(d, t);
        }

        /// @brief Returns the squared distance of p to edge ab, and the parametric coordinate of the closest point
        Eigen::Vector2d project_to_edge(
            const Eigen::Vector3d &a,
            const Eigen::Vector3d &b,
            const Eigen::Vector3d &p)
        {
            Eigen::Vector3d s = b - a;
            double t = (p - a).dot(s) / s.squaredNorm();
            t = std::min(1., std::max(t, 0.));
            double d = (p - (a + t * s)).squaredNorm();

            return Eigen::Vector2d(d, t);
        }

        void floydWarshall(const Eigen::MatrixXi &G, Eigen::MatrixXi &dist, Eigen::MatrixXi &parents)
        {
            int N = G.rows();
            dist = G;
            parents = -Eigen::MatrixXi::Ones(N, N);
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    if (dist(i, j) > 0 && dist(i, j) <= N)
                        parents(i, j) = i;

            for (int k = 0; k < N; k++)
            {
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        if (k == i || k == j || i == j)
                            continue;
                        if (dist(i, j) > dist(i, k) + dist(k, j))
                        {
                            dist(i, j) = dist(i, k) + dist(k, j);
                            parents(i, j) = parents(k, j);
                        }
                    }
                }
            }

            // validate
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (i == j)
                        continue;
                    const int p = parents(i, j);
                    if (dist(i, p) + G(p, j) != dist(i, j))
                        log_and_throw_error("[floydWarshall] Wrong closest distance!");

                    int cur = j;
                    while (parents(i, cur) != i)
                        cur = parents(i, cur);
                    if (parents(i, cur) != i)
                        log_and_throw_error("[floydWarshall] Wrong closest distance!");
                }
            }
        }

        bool are_same_edges(const Eigen::MatrixXi &A, const Eigen::MatrixXi &B)
        {
            if (A.rows() != B.rows())
                return false;

            for (int i = 0; i < A.rows(); i++)
            {
                bool flag = false;
                for (int j = 0; j < B.rows(); j++)
                {
                    if ((std::min(A(i, 0), A(i, 1)) == std::min(B(j, 0), B(j, 1)))
                    && (std::max(A(i, 0), A(i, 1)) == std::max(B(j, 0), B(j, 1))))
                    {
                        flag = true;
                        break;
                    }
                }
                if (!flag)
                    return false;
            }

            return true;
        }

        bool is_end_node(const Eigen::MatrixXi &edges, const int vid)
        {
            int cnt = 0;
            for (int i = 0; i < edges.size(); i++)
            {
                if (edges(i) == vid)
                    cnt++;
            }
            if (cnt == 0)
                log_and_throw_error("vid not found in is_end_node()!");
            return (cnt == 1);
        }

        void explode_trimesh(
            const Eigen::MatrixXd &Vin,
            const Eigen::MatrixXi &Fin,
            Eigen::MatrixXd &Vout,
            Eigen::MatrixXi &Fout,
            Eigen::VectorXi &index_map)
        {
            Vout = Eigen::MatrixXd::Zero(Fin.size(), 3);
            Fout = Fin;
            index_map.setZero(Vout.rows());
            for (int f = 0; f < Fin.rows(); f++)
            {
                for (int i = 0; i < Fin.cols(); i++)
                {
                    Vout.row(f * Fin.cols() + i) = Vin.row(Fin(f, i));
                    Fout(f, i) = f * Fin.cols() + i;
                    index_map(f * Fin.cols() + i) = Fin(f, i);
                }
            }
        }

    }

    void OBJMesh::read(const std::string &path)
    {
        igl::readOBJ(path, v, tc, cn, f, ftc, fn);
    }

    void OBJMesh::write(const std::string &path)
    {
        igl::writeOBJ(path, v, f, cn, fn, tc, ftc);
    }

    void GarmentSolver::save_result(
        const std::string &path,
        const int index,
        GarmentNLProblem &prob,
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::VectorXd &sol)
    {
        const Eigen::VectorXd full_disp = prob.reduced_to_full(sol);
        const Eigen::VectorXd complete_disp = prob.full_to_complete(full_disp);
        const Eigen::MatrixXd current_vertices = utils::unflatten(complete_disp, V.cols()) + V;

#ifdef POLYFEM_WITH_PARAVIEWO
        std::shared_ptr<paraviewo::ParaviewWriter> tmpw = std::make_shared<paraviewo::VTUWriter>();
        paraviewo::ParaviewWriter &writer = *tmpw;

        if (false)
        {
            Eigen::VectorXd total_grad = Eigen::VectorXd::Zero(complete_disp.size());
            std::unordered_set<std::string> existing_names;
            for (const auto &form : prob.forms())
            {
                Eigen::VectorXd grad;
                form->first_derivative(complete_disp, grad);
                std::string name = "grad_" + form->name();
                while (existing_names.count(name) != 0)
                    name += "_";
                existing_names.insert(name);
                grad.head(nc_avatar_v.rows() * 3).setZero();
                writer.add_field(name, utils::unflatten(grad, 3));
                total_grad += grad;
            }
            for (const auto &form : prob.full_forms())
            {
                Eigen::VectorXd grad_full, grad;
                form->first_derivative(full_disp, grad_full);
                std::string name = "grad_" + form->name();
                while (existing_names.count(name) != 0)
                    name += "_";
                existing_names.insert(name);
                grad.setZero(total_grad.size());
                grad.tail(grad_full.size() - 1) = grad_full.tail(grad_full.size() - 1);
                writer.add_field(name, utils::unflatten(grad, 3));
                total_grad += grad;
            }
            total_grad.head(nc_avatar_v.rows() * 3).setZero();
            writer.add_field("grad", utils::unflatten(total_grad, 3));

            Eigen::VectorXd body_ids = Eigen::VectorXd::Zero(V.rows());
            body_ids.head(nc_avatar_v.rows()).array() = 1;
            writer.add_field("body_ids", body_ids);

            logger().debug("Save VTU to {}", path + "/step_" + std::to_string(index) + ".vtu");
            writer.write_mesh(path + "/step_" + std::to_string(index) + ".vtu", current_vertices, F);
        }
#endif
        garment.v = current_vertices.bottomRows(garment.v.rows());

        // Save garment with group information if available
        if (!garment_objects.empty())
        {
            // Create OBJData structure with updated vertices
            OBJData garment_output_data;
            garment_output_data.V.resize(garment.v.rows());
            for (int i = 0; i < garment.v.rows(); ++i)
            {
                garment_output_data.V[i] = {garment.v(i, 0), garment.v(i, 1), garment.v(i, 2)};
            }

            // Convert faces back to vector format
            garment_output_data.F.resize(garment.f.rows());
            for (int i = 0; i < garment.f.rows(); ++i)
            {
                garment_output_data.F[i] = {garment.f(i, 0), garment.f(i, 1), garment.f(i, 2)};
            }

            // Copy group structure
            garment_output_data.objects = garment_objects;
            garment_output_data.face_to_group = garment_face_to_group;
            garment_output_data.face_to_object = garment_face_to_object;
            garment_output_data.mtl_filename = garment.mtl_filename;

            // Write with groups
            if (io::OBJWriter::write_with_groups(path + "/step_garment_" + std::to_string(index) + ".obj", garment_output_data))
            {
                logger().debug("Saved garment OBJ with groups: {}", path + "/step_garment_" + std::to_string(index) + ".obj");
            }
            else
            {
                logger().warn("Failed to save garment with groups, falling back to basic OBJ");
                garment.write(path + "/step_garment_" + std::to_string(index) + ".obj");
            }
        }
        else
        {
            garment.write(path + "/step_garment_" + std::to_string(index) + ".obj");
            logger().debug("Save OBJ to {}", path + "/step_garment_" + std::to_string(index) + ".obj");
        }

        // Write garment MTL file if materials exist
        if (!garment_materials.empty() && !garment.mtl_filename.empty())
        {
            std::string garment_mtl_dest = path + "/step_garment_" + std::to_string(index) + ".mtl";
            std::string source_dir = garment_mtl_source_path.empty() ? "" :
                                   std::filesystem::path(garment_mtl_source_path).parent_path().string();
            if (write_mtl_file(garment_mtl_dest, garment_materials, source_dir))
            {
                logger().debug("Saved garment MTL file: {}", garment_mtl_dest);
            }
        }

        // Save avatar with group information if available
        if (!avatar_objects.empty())
        {
            // Create OBJData structure with updated avatar vertices
            OBJData avatar_output_data;
            const Eigen::MatrixXd avatar_vertices = current_vertices.topRows(nc_avatar_v.rows());

            avatar_output_data.V.resize(avatar_vertices.rows());
            for (int i = 0; i < avatar_vertices.rows(); ++i)
            {
                avatar_output_data.V[i] = {avatar_vertices(i, 0), avatar_vertices(i, 1), avatar_vertices(i, 2)};
            }

            // Convert faces back to vector format
            avatar_output_data.F.resize(nc_avatar_f.rows());
            for (int i = 0; i < nc_avatar_f.rows(); ++i)
            {
                avatar_output_data.F[i] = {nc_avatar_f(i, 0), nc_avatar_f(i, 1), nc_avatar_f(i, 2)};
            }

            // Copy group structure
            avatar_output_data.objects = avatar_objects;
            avatar_output_data.face_to_group = avatar_face_to_group;
            avatar_output_data.face_to_object = avatar_face_to_object;
            avatar_output_data.mtl_filename = avatar_mtl_filename;

            // Copy texture coordinates (preserved from original avatar)
            avatar_output_data.VT = avatar_VT;
            avatar_output_data.FT = avatar_FT;

            // Skip writing normals to avoid artifacts from deformed geometry
            // Rendering engines will compute correct normals automatically from the mesh
            logger().debug("Skipping normals in output - will be computed automatically by rendering engines");

            // Write with groups
            if (io::OBJWriter::write_with_groups(path + "/step_avatar_" + std::to_string(index) + ".obj", avatar_output_data))
            {
                logger().debug("Saved avatar OBJ with groups: {}", path + "/step_avatar_" + std::to_string(index) + ".obj");
            }
            else
            {
                logger().warn("Failed to save avatar with groups, falling back to basic OBJ");
                igl::write_triangle_mesh(path + "/step_avatar_" + std::to_string(index) + ".obj", avatar_vertices, nc_avatar_f);
            }
        }
        else
        {
            igl::write_triangle_mesh(path + "/step_avatar_" + std::to_string(index) + ".obj", current_vertices.topRows(nc_avatar_v.rows()), nc_avatar_f);
        }

        // Write avatar MTL file if materials exist
        if (!avatar_materials.empty() && !avatar_mtl_filename.empty())
        {
            std::string avatar_mtl_dest = path + "/step_avatar_" + std::to_string(index) + ".mtl";
            std::string source_dir = avatar_mtl_source_path.empty() ? "" :
                                   std::filesystem::path(avatar_mtl_source_path).parent_path().string();
            if (write_mtl_file(avatar_mtl_dest, avatar_materials, source_dir))
            {
                logger().debug("Saved avatar MTL file: {}", avatar_mtl_dest);
            }
        }
    }

    Eigen::Vector3d bbox_size(const Eigen::Matrix<double, -1, 3> &V)
    {
        return V.colwise().maxCoeff() - V.colwise().minCoeff();
    }

    void GarmentSolver::load_garment_mesh(
        const std::string &mesh_path,
        const std::string &no_fit_spec_path)
	{
        // Read garment mesh with group information
        OBJData garment_obj_data;
        if (io::OBJReader::read_with_groups(mesh_path, garment_obj_data))
        {
            // Convert to Eigen matrices
            garment.v.resize(garment_obj_data.V.size(), 3);
            for (size_t i = 0; i < garment_obj_data.V.size(); ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    garment.v(i, j) = garment_obj_data.V[i][j];
                }
            }

            // Convert faces (assuming triangular faces)
            std::vector<std::vector<int>> triangular_faces;
            for (const auto &face : garment_obj_data.F)
            {
                if (face.size() == 3)
                {
                    triangular_faces.push_back(face);
                }
                else if (face.size() > 3)
                {
                    // Triangulate polygon faces
                    for (size_t i = 1; i < face.size() - 1; ++i)
                    {
                        triangular_faces.push_back({face[0], face[i], face[i + 1]});
                    }
                }
            }

            garment.f.resize(triangular_faces.size(), 3);
            for (size_t i = 0; i < triangular_faces.size(); ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    garment.f(i, j) = triangular_faces[i][j];
                }
            }

            // Store group information
            garment_objects = garment_obj_data.objects;
            garment_face_to_group = garment_obj_data.face_to_group;
            garment_face_to_object = garment_obj_data.face_to_object;
            garment.mtl_filename = garment_obj_data.mtl_filename;

            // Load MTL file if present
            if (!garment.mtl_filename.empty())
            {
                std::filesystem::path garment_dir = std::filesystem::path(mesh_path).parent_path();
                garment_mtl_source_path = (garment_dir / garment.mtl_filename).string();
                garment_materials = parse_mtl_file(garment_mtl_source_path);
                logger().info("Loaded garment with {} objects, {} groups, MTL file: {} with {} materials",
                             garment_objects.size(),
                             std::accumulate(garment_objects.begin(), garment_objects.end(), 0,
                                           [](int sum, const OBJObject& obj) { return sum + obj.groups.size(); }),
                             garment.mtl_filename, garment_materials.size());
            }
            else
            {
                logger().info("Loaded garment with {} objects, {} groups (no MTL file)",
                             garment_objects.size(),
                             std::accumulate(garment_objects.begin(), garment_objects.end(), 0,
                                           [](int sum, const OBJObject& obj) { return sum + obj.groups.size(); }));
            }
        }
        else
        {
            // Fallback to basic reading if group reading fails
            logger().warn("Failed to read garment with groups, falling back to basic OBJ reading");
            garment.read(mesh_path);
        }

        if (std::filesystem::exists(no_fit_spec_path))
        {
            Eigen::MatrixXi tmp_vids;
            io::read_matrix<int>(no_fit_spec_path, tmp_vids);

            if (tmp_vids.maxCoeff() >= garment.v.rows() || tmp_vids.minCoeff() < 0)
                log_and_throw_error("Vertex ID {} in no-fit.txt out of range!");

            Eigen::VectorXi vmask = Eigen::VectorXi::Zero(garment.v.rows());
            for (int i = 0; i < tmp_vids.size(); i++)
                vmask(tmp_vids(i)) = 1;

            for (int i = 0; i < garment.f.rows(); i++)
                if (vmask(garment.f(i, 0)) && vmask(garment.f(i, 1)) && vmask(garment.f(i, 2)))
                    not_fit_fids.push_back(i);
        }
        else
            logger().debug("Cannot find {}, will fit the garment tightly everywhere...", no_fit_spec_path);

        // if (std::filesystem::exists(path + "/skin.txt"))
        // {
        //     log_and_throw_error("Utilizing garment skinning weight is not supported!");

        //     io::read_matrix(path + "/skin.txt", garment_skinning_weights);
        //     assert(garment_skinning_weights.rows() == skeleton_v.rows());
        //     assert(garment.v.rows() == garment_skinning_weights.cols());
        //     assert(garment_skinning_weights.minCoeff() >= 0. && garment_skinning_weights.maxCoeff() <= 1.);
        // }
        // else
        // if (n_refs > 0) {
        //     while (n_refs-- > 0)
        //     {
        //         std::tie(garment.v, garment.f) = refine(garment.v, garment.f);

        //         std::vector<int> not_fit_fids_new;
        //         for (int i = 0; i < not_fit_fids.size(); i++)
        //             for (int j = 0; j < 4; j++)
        //                 not_fit_fids_new.push_back(not_fit_fids[i] * 4 + j);
        //         std::swap(not_fit_fids, not_fit_fids_new);
        //     }
        //     assert(n_refs == 0);
        // }

		// remove duplicate vertices in the garment
        // remove_duplicate_vertices(garment.v, garment.f, 1e-6);
	}

    void GarmentSolver::check_intersections(
        const ipc::CollisionMesh &collision_mesh,
        const Eigen::MatrixXd &collision_vertices) const
    {
        auto ids = ipc::my_has_intersections(collision_mesh, collision_vertices, ipc::BroadPhaseMethod::BVH);
        if (ids[0] >= 0)
        {
            // Create basic OBJData for intersection debugging
            OBJData intersection_data;
            intersection_data.V.resize(collision_vertices.rows());
            for (int i = 0; i < collision_vertices.rows(); ++i)
            {
                intersection_data.V[i] = {collision_vertices(i, 0), collision_vertices(i, 1), collision_vertices(i, 2)};
            }

            // Add faces
            intersection_data.F.resize(collision_mesh.faces().rows());
            for (int i = 0; i < collision_mesh.faces().rows(); ++i)
            {
                intersection_data.F[i] = {collision_mesh.faces()(i, 0), collision_mesh.faces()(i, 1), collision_mesh.faces()(i, 2)};
            }

            // Create default object and group
            OBJObject default_obj;
            default_obj.name = "intersection";
            OBJGroup default_group;
            default_group.name = "default";
            for (int i = 0; i < intersection_data.F.size(); ++i)
            {
                default_group.face_indices.push_back(i);
            }
            default_obj.groups.push_back(default_group);
            intersection_data.objects.push_back(default_obj);

            // Set face mappings
            intersection_data.face_to_object.resize(intersection_data.F.size(), 0);
            intersection_data.face_to_group.resize(intersection_data.F.size(), 0);

            io::OBJWriter::write_with_groups(out_folder + "/intersection.obj", intersection_data);

            // Create intersecting pair data
            OBJData pair_data;
            pair_data.V = intersection_data.V; // Same vertices

            // Add the specific intersecting edge and face
            pair_data.F.push_back({ids[2], ids[3], ids[4]}); // Face

            OBJObject pair_obj;
            pair_obj.name = "intersecting_pair";
            OBJGroup pair_group;
            pair_group.name = "default";
            pair_group.face_indices.push_back(0);
            pair_obj.groups.push_back(pair_group);
            pair_data.objects.push_back(pair_obj);

            pair_data.face_to_object.push_back(0);
            pair_data.face_to_group.push_back(0);

            io::OBJWriter::write_with_groups(out_folder + "/intersecting_pair.obj", pair_data);
            log_and_throw_error("Unable to solve, initial solution has intersections!");
        }
    }

    void GarmentSolver::read_meshes(
        const std::string &avatar_mesh_path,
        const std::string &source_skeleton_path,
        const std::string &target_skeleton_path,
        const std::string &target_avatar_skinning_weights_path)
    {
        // Read avatar mesh with group information
        OBJData avatar_obj_data;
        if (io::OBJReader::read_with_groups(avatar_mesh_path, avatar_obj_data))
        {
            // Convert to Eigen matrices
            avatar_v.resize(avatar_obj_data.V.size(), 3);
            for (size_t i = 0; i < avatar_obj_data.V.size(); ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    avatar_v(i, j) = avatar_obj_data.V[i][j];
                }
            }

            // Convert faces (assuming triangular faces)
            std::vector<std::vector<int>> triangular_faces;
            for (const auto &face : avatar_obj_data.F)
            {
                if (face.size() == 3)
                {
                    triangular_faces.push_back(face);
                }
                else if (face.size() > 3)
                {
                    // Triangulate polygon faces
                    for (size_t i = 1; i < face.size() - 1; ++i)
                    {
                        triangular_faces.push_back({face[0], face[i], face[i + 1]});
                    }
                }
            }

            avatar_f.resize(triangular_faces.size(), 3);
            for (size_t i = 0; i < triangular_faces.size(); ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    avatar_f(i, j) = triangular_faces[i][j];
                }
            }

            // Store group information
            avatar_objects = avatar_obj_data.objects;
            avatar_face_to_group = avatar_obj_data.face_to_group;
            avatar_face_to_object = avatar_obj_data.face_to_object;
            avatar_mtl_filename = avatar_obj_data.mtl_filename;

            // Store texture coordinates and normals
            avatar_VT = avatar_obj_data.VT;
            avatar_VN = avatar_obj_data.VN;
            avatar_FT = avatar_obj_data.FT;
            avatar_FN = avatar_obj_data.FN;

            // Load MTL file if present
            if (!avatar_mtl_filename.empty())
            {
                std::filesystem::path avatar_dir = std::filesystem::path(avatar_mesh_path).parent_path();
                avatar_mtl_source_path = (avatar_dir / avatar_mtl_filename).string();
                avatar_materials = parse_mtl_file(avatar_mtl_source_path);
                logger().info("Loaded avatar with {} objects, {} groups, MTL file: {} with {} materials",
                             avatar_objects.size(),
                             std::accumulate(avatar_objects.begin(), avatar_objects.end(), 0,
                                           [](int sum, const OBJObject& obj) { return sum + obj.groups.size(); }),
                             avatar_mtl_filename, avatar_materials.size());
            }
            else
            {
                logger().info("Loaded avatar with {} objects, {} groups (no MTL file)",
                             avatar_objects.size(),
                             std::accumulate(avatar_objects.begin(), avatar_objects.end(), 0,
                                           [](int sum, const OBJObject& obj) { return sum + obj.groups.size(); }));
            }
        }
        else
        {
            // Fallback to basic reading if group reading fails
            logger().warn("Failed to read avatar with groups, falling back to basic OBJ reading");
            igl::read_triangle_mesh(avatar_mesh_path, avatar_v, avatar_f);
        }

        read_edge_mesh(source_skeleton_path, skeleton_v, skeleton_b);
        read_edge_mesh(target_skeleton_path, target_skeleton_v, target_skeleton_b);
        if (!are_same_edges(skeleton_b, target_skeleton_b))
            log_and_throw_error("Inconsistent skeletons!");
        target_skeleton_b = skeleton_b;

        if (std::filesystem::exists(target_avatar_skinning_weights_path))
        {
            io::read_matrix(target_avatar_skinning_weights_path, target_avatar_skinning_weights);
            if (target_avatar_skinning_weights.rows() != skeleton_v.rows()
                || avatar_v.rows() != target_avatar_skinning_weights.cols())
                log_and_throw_error("Inconsistent skin weights dimension with the number of vertices and bones! Skin weights dimension: {}x{}, number of bones: {}, number of vertices: {}", target_avatar_skinning_weights.rows(), target_avatar_skinning_weights.cols(), skeleton_v.rows(), avatar_v.rows());
        }
        else
        {
            target_avatar_skinning_weights.setZero(0, 0);
            logger().warn("Cannot find target avatar skinning weights, use pure distance-based projection instead...");
        }
    }

    void GarmentSolver::project_avatar_to_skeleton()
    {
        Eigen::MatrixXi graph, shared_vtx, dist, parent;
        {
            graph = Eigen::MatrixXi::Ones(skeleton_b.rows(), skeleton_b.rows()) * (skeleton_b.rows() + 1);
            shared_vtx = -Eigen::MatrixXi::Ones(skeleton_b.rows(), skeleton_b.rows());
            for (int i = 0; i < skeleton_b.rows(); i++)
            {
                graph(i, i) = 0;
                for (int j = 0; j < skeleton_b.rows(); j++)
                {
                    bool adjacent = (skeleton_b(i, 0) == skeleton_b(j, 0)) ||
                                    (skeleton_b(i, 0) == skeleton_b(j, 1)) ||
                                    (skeleton_b(i, 1) == skeleton_b(j, 0)) ||
                                    (skeleton_b(i, 1) == skeleton_b(j, 1));
                    if (i != j && adjacent)
                    {
                        graph(i, j) = 1;
                        if ((skeleton_b(i, 0) == skeleton_b(j, 0)) || (skeleton_b(i, 0) == skeleton_b(j, 1)))
                            shared_vtx(i, j) = skeleton_b(i, 0);
                        else
                            shared_vtx(i, j) = skeleton_b(i, 1);
                    }
                }
            }
            floydWarshall(graph, dist, parent);
        }

        const bool has_target_avatar_skin_weights = target_avatar_skinning_weights.size() > 0;

        // Keep original mesh connectivity instead of exploding
        nc_avatar_v = avatar_v;
        nc_avatar_f = avatar_f;

        // IMPORTANT: Preserve texture coordinates during projection, but recompute normals
        // The avatar projection dramatically changes vertex positions, making original normals incorrect
        logger().info("Preserving texture coordinates and material data during avatar projection");
        logger().info("Normals will be recomputed after projection to match the new geometry");

        Eigen::MatrixXd new_skinning_weights;
        if (has_target_avatar_skin_weights)
            new_skinning_weights = target_avatar_skinning_weights;


        Eigen::VectorXi eid;
        Eigen::VectorXd coord;
        Eigen::MatrixXd skinny_avatar_v_debug;
        // first pass
        {
            const int N = nc_avatar_v.rows();
            Eigen::VectorXd dists(N);
            dists.setConstant(std::numeric_limits<double>::max());
            coord.setZero(N);
            eid = -Eigen::VectorXi::Ones(N);
            for (int i = 0; i < N; i++)
            {
                Eigen::Index maxRow = 0;
                if (has_target_avatar_skin_weights)
                {
                    Eigen::Index maxCol;
                    const double max_skin_weight = new_skinning_weights.col(i).maxCoeff(&maxRow, &maxCol);
                    assert(maxCol == 0);
                }

                for (int e = 0; e < target_skeleton_b.rows(); e++)
                {
                    if (!has_target_avatar_skin_weights || target_skeleton_b(e, 0) == maxRow || target_skeleton_b(e, 1) == maxRow)
                    {
                        Eigen::Vector2d tmp1 = project_to_edge(target_skeleton_v.row(target_skeleton_b(e, 0)), target_skeleton_v.row(target_skeleton_b(e, 1)), nc_avatar_v.row(i));
                        if (tmp1(0) < dists(i))
                        {
                            dists(i) = tmp1(0);
                            {
                                Eigen::Vector2d tmp2 = project_to_line(target_skeleton_v.row(target_skeleton_b(e, 0)), target_skeleton_v.row(target_skeleton_b(e, 1)), nc_avatar_v.row(i));
                                if (!is_end_node(target_skeleton_b, target_skeleton_b(e, 0)))
                                    tmp2(1) = std::max(0., tmp2(1));
                                if (!is_end_node(target_skeleton_b, target_skeleton_b(e, 1)))
                                    tmp2(1) = std::min(1., tmp2(1));
                                coord(i) = tmp2(1);
                            }
                            eid(i) = e;
                        }
                    }
                }
                if (eid(i) < 0)
                    log_and_throw_error("Failed to project vertex to the bone!");
            }

            skinny_avatar_v.setZero(nc_avatar_v.rows(), nc_avatar_v.cols());
            for (int i = 0; i < nc_avatar_v.rows(); i++)
                skinny_avatar_v.row(i) += coord(i) * (skeleton_v(skeleton_b(eid(i), 1), Eigen::all) - skeleton_v(skeleton_b(eid(i), 0), Eigen::all)) + skeleton_v(skeleton_b(eid(i), 0), Eigen::all);

            skinny_avatar_v_debug.setZero(nc_avatar_v.rows(), nc_avatar_v.cols());
            for (int i = 0; i < nc_avatar_v.rows(); i++)
                skinny_avatar_v_debug.row(i) += coord(i) * (target_skeleton_v(skeleton_b(eid(i), 1), Eigen::all) - target_skeleton_v(skeleton_b(eid(i), 0), Eigen::all)) + target_skeleton_v(skeleton_b(eid(i), 0), Eigen::all);

            skinny_avatar_f = nc_avatar_f;
        }

        // igl::write_triangle_mesh(out_folder + "/avatar_old.obj", nc_avatar_v, nc_avatar_f);
        // igl::write_triangle_mesh(out_folder + "/projected_avatar_old_source.obj", skinny_avatar_v, nc_avatar_f);
        // igl::write_triangle_mesh(out_folder + "/projected_avatar_old_target.obj", skinny_avatar_v_debug, nc_avatar_f);

        // Skip the iterative distance reduction that assumes exploded mesh
        // This preserves the original mesh topology
        logger().info("Preserving mesh topology - skipping iterative distance reduction");

        {
            Eigen::MatrixXd tmp_v(nc_avatar_v.rows(), nc_avatar_v.cols() + skinny_avatar_v.cols());
            tmp_v << nc_avatar_v, skinny_avatar_v;
            const auto [svi, svj] = remove_duplicate_vertices(tmp_v, nc_avatar_f, 1e-10);

            nc_avatar_v = tmp_v.template leftCols<3>();
            skinny_avatar_v = tmp_v.template rightCols<3>();
            skinny_avatar_f = nc_avatar_f;
        }

        skinny_avatar_v += (nc_avatar_v - skinny_avatar_v) * 1e-2;
    }

    void GarmentSolver::normalize_meshes()
    {
        // Center offset
        const Eigen::RowVector3d center_offset = skeleton_v.colwise().sum() / skeleton_v.rows();
        skeleton_v.rowwise() -= center_offset;
        garment.v.rowwise() -= center_offset;

        // Source side
        const double source_scaling = 2. / bbox_size(skeleton_v).maxCoeff();
        skeleton_v *= source_scaling;
        garment.v *= source_scaling;
        // skinny_avatar_v *= source_scaling;

        // Target side
        const double target_scaling = bbox_size(skeleton_v).maxCoeff() / bbox_size(target_skeleton_v).maxCoeff();
        const Eigen::Vector3d center = skeleton_v.colwise().sum() / skeleton_v.rows() - target_scaling * avatar_v.colwise().sum() / avatar_v.rows();
        Transformation<3> trans(target_scaling * Eigen::Matrix3d::Identity(), center);

        trans.apply(avatar_v);
        trans.apply(target_skeleton_v);
    }

	json init(const json &p_args_in, const bool strict_validation)
	{
		json args_in = p_args_in; // mutable copy
        json args;

		utils::apply_common_params(args_in);

		// CHECK validity json
		json rules;
		jse::JSE jse;
		{
			jse.strict = strict_validation;
			const std::string polyfem_input_spec = POLYFEM_INPUT_SPEC;
			std::ifstream file(polyfem_input_spec);

			if (file.is_open())
				file >> rules;
			else
			{
				logger().error("unable to open {} rules", polyfem_input_spec);
				throw std::runtime_error("Invald spec file");
			}

			jse.include_directories.push_back(POLYFEM_JSON_SPEC_DIR);
			jse.include_directories.push_back(POLYSOLVE_JSON_SPEC_DIR);
			rules = jse.inject_include(rules);

			polysolve::linear::Solver::apply_default_solver(rules, "/solver/linear");
		}

		polysolve::linear::Solver::select_valid_solver(args_in["solver"]["linear"], logger());

		// Use the /solver/nonlinear settings as the default for /solver/augmented_lagrangian/nonlinear
		if (args_in.contains("/solver/nonlinear"_json_pointer))
		{
			if (args_in.contains("/solver/augmented_lagrangian/nonlinear"_json_pointer))
			{
				assert(args_in["solver"]["augmented_lagrangian"]["nonlinear"].is_object());
				// Merge the augmented lagrangian settings into the nonlinear settings,
				// and then replace the augmented lagrangian settings with the merged settings.
				json nonlinear = args_in["solver"]["nonlinear"]; // copy
				nonlinear.merge_patch(args_in["solver"]["augmented_lagrangian"]["nonlinear"]);
				args_in["solver"]["augmented_lagrangian"]["nonlinear"] = nonlinear;
			}
			else
			{
				// Copy the nonlinear settings to the augmented_lagrangian settings
				args_in["solver"]["augmented_lagrangian"]["nonlinear"] = args_in["solver"]["nonlinear"];
			}
		}

		const bool valid_input = jse.verify_json(args_in, rules);

		if (!valid_input)
		{
			logger().error("invalid input json:\n{}", jse.log2str());
			throw std::runtime_error("Invald input json file");
		}
		// end of check

		args = jse.inject_defaults(args_in, rules);

		// Save output directory and resolve output paths dynamically
		const std::string output_dir = utils::resolve_path(args["output"]["directory"],
            utils::is_param_valid(args, "root_path") ? args["root_path"].get<std::string>() : "",
            false);

		if (!output_dir.empty())
		{
			std::filesystem::create_directories(output_dir);
		}

        // set logger
        {
            spdlog::level::level_enum log_level = args["output"]["log"]["level"];

            spdlog::set_level(log_level);
			logger().set_level(log_level);
			ipc::logger().set_level(log_level);

            spdlog::flush_every(std::chrono::seconds(3));
        }

		logger().info("Saving output to {}", output_dir);

		const unsigned int thread_in = args["solver"]["max_threads"];
		utils::NThread::get().set_num_threads(thread_in);

        return args;
	}
}
