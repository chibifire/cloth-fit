#include "optimize.hpp"

#include <polyfem/solver/GarmentNLProblem.hpp>
#include <polyfem/solver/forms/ContactForm.hpp>
#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/mesh/MeshUtils.hpp>

#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/write_triangle_mesh.h>

#include <paraviewo/ParaviewWriter.hpp>
#include <paraviewo/VTUWriter.hpp>

#include <ipc/ipc.hpp>
#include <ipc/distance/point_edge.hpp>

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

        // explode avatar mesh
        Eigen::MatrixXd new_avatar_v(avatar_f.size(), 3);
        for (int f = 0; f < avatar_f.rows(); f++)
        {
            for (int i = 0; i < avatar_f.cols(); i++)
            {
                new_avatar_v.row(f * avatar_f.cols() + i) = avatar_v.row(avatar_f(f, i));
                avatar_f(f, i) = f * avatar_f.cols() + i;
            }
        }
        std::swap(new_avatar_v, avatar_v);

        Eigen::VectorXi eid;
        Eigen::VectorXd coord;
        // first pass
        {
            std::tie(eid, coord) = project_to_edge_mesh(target_skeleton_v, target_skeleton_b, avatar_v);

            skinny_avatar_v.setZero(avatar_v.rows(), avatar_v.cols());
            for (int i = 0; i < avatar_v.rows(); i++)
                skinny_avatar_v.row(i) += coord(i) * (skeleton_v(skeleton_b(eid(i), 1), Eigen::all) - skeleton_v(skeleton_b(eid(i), 0), Eigen::all)) + skeleton_v(skeleton_b(eid(i), 0), Eigen::all);
            skinny_avatar_f = avatar_f;
        }

        igl::write_triangle_mesh(out_folder + "/projected_avatar_old.obj", skinny_avatar_v, avatar_f);

        // iteratively reduce distance
        for (int iter = 0; iter < 10; iter++) {
            const int n_faces = skinny_avatar_f.rows();
            std::vector<Eigen::Matrix<double, 6, 3>> new_faces;
            for (int f = 0; f < n_faces; f++)
            {
                int max_dist = 0;
                int max_dist_le = -1;
                for (int le = 0; le < 3; le++)
                {
                    const int a = f * 3 + le;
                    const int b = f * 3 + (le + 1) % 3;
                    const int c = f * 3 + (le + 2) % 3;

                    bool edge_overlap_with_skeleton = false;
                    for (int i = 0; i < target_skeleton_b.rows(); i++) 
                    {
                        Eigen::Vector3d vb = skinny_avatar_v.row(b);
                        Eigen::Vector3d va = skinny_avatar_v.row(a);
                        Eigen::Vector3d b0 = skeleton_v.row(target_skeleton_b(i, 0));
                        Eigen::Vector3d b1 = skeleton_v.row(target_skeleton_b(i, 1));

                        if (ipc::point_edge_distance(vb, b0, b1) < 1e-4 * (b1 - b0).squaredNorm() &&
                            ipc::point_edge_distance(va, b0, b1) < 1e-4 * (b1 - b0).squaredNorm())
                        {
                            edge_overlap_with_skeleton = true;
                            break;
                        }
                    }

                    int source = eid(b);
                    int cur = eid(a);

                    if (edge_overlap_with_skeleton || source == cur)
                        continue;
                    
                    if (dist(cur, source) > max_dist)
                    {
                        max_dist = dist(cur, source);
                        max_dist_le = le;
                    }
                }

                for (int le = 0; le < 3; le++)
                {
                    if (max_dist_le != le)
                        continue;
                    
                    const int a = f * 3 + le;
                    const int b = f * 3 + (le + 1) % 3;
                    const int c = f * 3 + (le + 2) % 3;

                    int source = eid(b);
                    int cur = eid(a);

                    std::vector<std::array<int, 2>> inserted_tmp;
                    while (cur != source)
                    {
                        if (cur < 0 || source < 0 || source >= parent.rows() || cur >= parent.cols() || parent(source, cur) < 0)
                            std::cout << std::endl;
                        std::array<int, 2> tmp{shared_vtx(cur, parent(source, cur)), cur};
                        inserted_tmp.push_back(tmp);
                        cur = parent(source, cur);
                    }

                    Eigen::RowVector3d p_prev = skinny_avatar_v.row(a);
                    for (int k = 0; k < inserted_tmp.size(); k++)
                    {
                        Eigen::RowVector3d p = skeleton_v.row(inserted_tmp[k][0]);

                        Eigen::Matrix<double, 6, 3> X;
                        X << p_prev, p, skinny_avatar_v.row(c),
                             avatar_v.row(a) + (avatar_v.row(b) - avatar_v.row(a)) * ((double)k / (inserted_tmp.size() + 1)),
                             avatar_v.row(a) + (avatar_v.row(b) - avatar_v.row(a)) * ((double)(k+1) / (inserted_tmp.size() + 1)),
                             avatar_v.row(c);
                        new_faces.push_back(X);
                        eid.conservativeResize(eid.size() + 3);
                        eid.tail(3) << (int)((k == 0) ? eid(a) : inserted_tmp[k-1][1]), inserted_tmp[k][1], eid(c);

                        p_prev = p;
                    }

                    skinny_avatar_v.row(a) = p_prev;
                    avatar_v.row(a) += (avatar_v.row(b) - avatar_v.row(a)) * ((double)inserted_tmp.size() / (inserted_tmp.size() + 1));
                    eid(a) = inserted_tmp.back()[1];
                }
            }

            if (new_faces.size() == 0)
                break;
        
            {
                Eigen::MatrixXd tmp(skinny_avatar_v.rows() + new_faces.size() * 3, 3);
                tmp.topRows(skinny_avatar_v.rows()) = skinny_avatar_v;
                for (int i = 0; i < new_faces.size(); i++)
                    tmp.block(skinny_avatar_v.rows() + 3 * i, 0, 3, 3) = new_faces[i].topRows(3);
                std::swap(skinny_avatar_v, tmp);
            }

            {
                Eigen::MatrixXd tmp(avatar_v.rows() + new_faces.size() * 3, 3);
                tmp.topRows(avatar_v.rows()) = avatar_v;
                for (int i = 0; i < new_faces.size(); i++)
                    tmp.block(avatar_v.rows() + 3 * i, 0, 3, 3) = new_faces[i].bottomRows(3);
                std::swap(avatar_v, tmp);
            }

            skinny_avatar_f = Eigen::VectorXi::LinSpaced(skinny_avatar_v.rows(), 0, skinny_avatar_v.rows() - 1).reshaped(3, skinny_avatar_v.rows() / 3).transpose();
            avatar_f = skinny_avatar_f;
        
            igl::write_triangle_mesh(out_folder + "/projected_avatar_new_" + std::to_string(iter) + ".obj", skinny_avatar_v, skinny_avatar_f);
            igl::write_triangle_mesh(out_folder + "/avatar_new_" + std::to_string(iter) + ".obj", avatar_v, avatar_f);
        }

        {
            Eigen::MatrixXi sf;
            Eigen::MatrixXd sv;
            Eigen::VectorXi svi, svj;
            igl::remove_duplicate_vertices(avatar_v, avatar_f, 1e-10, sv, svi, svj, sf);
            for (int i = 0; i < sf.rows(); i++)
            {
                if (sf(i, 0) == sf(i, 1) || sf(i, 2) == sf(i, 1) || sf(i, 0) == sf(i, 2))
                    log_and_throw_error("Treshold in igl::remove_duplicate_vertices is too large!!");
            }

            std::swap(sv, avatar_v);
            std::swap(sf, avatar_f);

            skinny_avatar_v = skinny_avatar_v(svi, Eigen::all).eval();
            skinny_avatar_f = avatar_f;
        }
        
        skinny_avatar_v += (avatar_v - skinny_avatar_v) * 1e-3;
    }

    void GarmentSolver::normalize_meshes()
    {
        // Source side
        const double source_scaling = 1e2;
        skeleton_v *= source_scaling;
        garment_v *= source_scaling;
        // skinny_avatar_v *= source_scaling;

        // Target side
        const double target_scaling = bbox_size(skeleton_v).maxCoeff() / bbox_size(target_skeleton_v).maxCoeff();
        const Eigen::Vector3d center = skeleton_v.colwise().sum() / skeleton_v.rows() - target_scaling * avatar_v.colwise().sum() / avatar_v.rows();
        Transformation<3> trans(target_scaling * Eigen::Matrix3d::Identity(), center);

        trans.apply(avatar_v);
        trans.apply(target_skeleton_v);
    }
}