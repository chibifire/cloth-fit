#include "optimize.hpp"

#include <polyfem/solver/GarmentNLProblem.hpp>
#include <polyfem/solver/forms/ContactForm.hpp>
#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/mesh/MeshUtils.hpp>

#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>

#include <paraviewo/ParaviewWriter.hpp>
#include <paraviewo/VTUWriter.hpp>

#include <ipc/ipc.hpp>

#include <unordered_set>

using namespace polyfem::solver;
using namespace polyfem::mesh;

namespace polyfem {
    namespace {
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

        std::tuple<Eigen::VectorXi, Eigen::VectorXd> 
        project_to_edge_mesh(
            const Eigen::MatrixXd &V,
            const Eigen::MatrixXi &E,
            const Eigen::MatrixXd &P)
        {
            const int N = P.rows();
            Eigen::VectorXd d(N), t(N);
            d.setConstant(std::numeric_limits<double>::max());
            t.setZero();
            Eigen::VectorXi I = -Eigen::VectorXi::Ones(N);
            for (int e = 0; e < E.rows(); e++)
            {
                for (int i = 0; i < N; i++)
                {
                    Eigen::Vector2d tmp = project_to_edge(V.row(E(e, 0)), V.row(E(e, 1)), P.row(i));
                    if (tmp(0) < d(i))
                    {
                        t(i) = tmp(1);
                        d(i) = tmp(0);
                        I(i) = e;
                    }
                }
            }

            return {I, t};
        }
    }
    void save_vtu(
        const std::string &path,
        GarmentNLProblem &prob,
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const int n_avatar_vertices,
        const Eigen::VectorXd &sol)
    {
        std::shared_ptr<paraviewo::ParaviewWriter> tmpw = std::make_shared<paraviewo::VTUWriter>();
        paraviewo::ParaviewWriter &writer = *tmpw;

        const Eigen::VectorXd complete_disp = prob.full_to_complete(prob.reduced_to_full(sol));

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
            writer.add_field(name, utils::unflatten(grad, 3));
            total_grad += grad;
        }
        writer.add_field("grad", utils::unflatten(total_grad, 3));

        Eigen::VectorXd body_ids = Eigen::VectorXd::Zero(V.rows());
        body_ids.head(n_avatar_vertices).array() = 1;
        writer.add_field("body_ids", body_ids);

        writer.write_mesh(path, utils::unflatten(complete_disp, V.cols()) + V, F);
    }

    Eigen::Vector3d bbox_size(const Eigen::Matrix<double, -1, 3> &V)
    {
        return V.colwise().maxCoeff() - V.colwise().minCoeff();
    }

    void GarmentSolver::load_garment_mesh(const std::string &path, int n_refs)
	{
		igl::read_triangle_mesh(path, garment_v, garment_f);

		// remove duplicate vertices in the garment
		{
			Eigen::VectorXi svi, svj;
			Eigen::MatrixXi sf;
			Eigen::MatrixXd sv;
			igl::remove_duplicate_vertices(garment_v, garment_f, 1e-4, sv, svi, svj, sf);
			std::swap(sv, garment_v);
			std::swap(sf, garment_f);
		}

		while (n_refs-- > 0)
			std::tie(garment_v, garment_f) = refine(garment_v, garment_f);
	}

    void GarmentSolver::check_intersections(
        const ipc::CollisionMesh &collision_mesh,
        const Eigen::MatrixXd &collision_vertices) const
    {
        auto ids = ipc::my_has_intersections(collision_mesh, collision_vertices, ipc::BroadPhaseMethod::BVH);
        if (ids[0] >= 0)
        {
            io::OBJWriter::write(
                out_folder + "/intersection.obj", collision_vertices,
                collision_mesh.edges(), collision_mesh.faces());
            Eigen::MatrixXi edge(1, 2);
            edge << ids[0], ids[1];
            Eigen::MatrixXi face(1, 3);
            face << ids[2], ids[3], ids[4];
            io::OBJWriter::write(
                out_folder + "/intersecting_pair.obj", collision_vertices,
                edge, face);
            log_and_throw_error("Unable to solve, initial solution has intersections!");
        }
    }
        
    void GarmentSolver::read_meshes(
        const std::string &avatar_mesh_path,
        const std::string &source_skeleton_path,
        const std::string &target_skeleton_path)
    {
        igl::read_triangle_mesh(avatar_mesh_path, avatar_v, avatar_f);
        read_edge_mesh(source_skeleton_path, skeleton_v, skeleton_b);
        read_edge_mesh(target_skeleton_path, target_skeleton_v, target_skeleton_b);
        assert((skeleton_b - target_skeleton_b).squaredNorm() < 1);
    }

    void GarmentSolver::project_avatar_to_skeleton()
    {
        const auto [eid, coord] = project_to_edge_mesh(target_skeleton_v, target_skeleton_b, avatar_v);

        // Eigen::MatrixXi E;
        // {
        //     skinny_avatar_v.setZero(avatar_v.rows(), avatar_v.cols());
        //     for (int i = 0; i < avatar_v.rows(); i++)
        //         skinny_avatar_v.row(i) += coord(i) * (target_skeleton_v(skeleton_b(eid(i), 1), Eigen::all) - target_skeleton_v(skeleton_b(eid(i), 0), Eigen::all)) + target_skeleton_v(skeleton_b(eid(i), 0), Eigen::all);

        //     io::OBJWriter::write("debug_project_A.obj", skinny_avatar_v, E, avatar_f);
        //     io::OBJWriter::write("debug_skeleton_A.obj", target_skeleton_v, target_skeleton_b, E);
        // }

        skinny_avatar_v.setZero(avatar_v.rows(), avatar_v.cols());
        for (int i = 0; i < avatar_v.rows(); i++)
            skinny_avatar_v.row(i) += coord(i) * (skeleton_v(skeleton_b(eid(i), 1), Eigen::all) - skeleton_v(skeleton_b(eid(i), 0), Eigen::all)) + skeleton_v(skeleton_b(eid(i), 0), Eigen::all);

        // io::OBJWriter::write("debug_project_B.obj", skinny_avatar_v, E, avatar_f);
        // io::OBJWriter::write("debug_skeleton_B.obj", skeleton_v, skeleton_b, E);

        skinny_avatar_v += (avatar_v - skinny_avatar_v) * 1e-4;
    }

    void GarmentSolver::normalize_meshes()
    {
        // Source side
        const double source_scaling = 1e2;
        skeleton_v *= source_scaling;
        garment_v *= source_scaling;
        skinny_avatar_v *= source_scaling;

        // Target side
        const double target_scaling = bbox_size(skinny_avatar_v).maxCoeff() / bbox_size(target_skeleton_v).maxCoeff();
        const Eigen::Vector3d center = (skinny_avatar_v.colwise().sum() - target_scaling * avatar_v.colwise().sum()) / avatar_v.rows();
        Transformation<3> trans(target_scaling * Eigen::Matrix3d::Identity(), center);

        trans.apply(avatar_v);
        trans.apply(target_skeleton_v);
    }
}