#include "optimize.hpp"

#include <polyfem/solver/GarmentNLProblem.hpp>
#include <polyfem/solver/forms/ContactForm.hpp>
#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/io/MatrixIO.hpp>
#include <polyfem/mesh/MeshUtils.hpp>

#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readOBJ.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/write_triangle_mesh.h>
#include <igl/writeOBJ.h>

#include <paraviewo/ParaviewWriter.hpp>
#include <paraviewo/VTUWriter.hpp>

#include <ipc/ipc.hpp>
#include <ipc/distance/point_edge.hpp>

#include <unordered_set>
#include <filesystem>

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
        std::shared_ptr<paraviewo::ParaviewWriter> tmpw = std::make_shared<paraviewo::VTUWriter>();
        paraviewo::ParaviewWriter &writer = *tmpw;

        const Eigen::VectorXd full_disp = prob.reduced_to_full(sol);
        const Eigen::VectorXd complete_disp = prob.full_to_complete(full_disp);

        // Eigen::VectorXd total_grad = Eigen::VectorXd::Zero(complete_disp.size());
        // std::unordered_set<std::string> existing_names;
        // for (const auto &form : prob.forms())
        // {
        //     Eigen::VectorXd grad;
        //     form->first_derivative(complete_disp, grad);
        //     std::string name = "grad_" + form->name();
        //     while (existing_names.count(name) != 0)
        //         name += "_";
        //     existing_names.insert(name);
        //     grad.head(n_avatar_vertices() * 3).setZero();
        //     writer.add_field(name, utils::unflatten(grad, 3));
        //     total_grad += grad;
        // }
        // for (const auto &form : prob.full_forms())
        // {
        //     Eigen::VectorXd grad_full, grad;
        //     form->first_derivative(full_disp, grad_full);
        //     std::string name = "grad_" + form->name();
        //     while (existing_names.count(name) != 0)
        //         name += "_";
        //     existing_names.insert(name);
        //     grad.setZero(total_grad.size());
        //     grad.tail(grad_full.size() - 1) = grad_full.tail(grad_full.size() - 1);
        //     writer.add_field(name, utils::unflatten(grad, 3));
        //     total_grad += grad;
        // }
        // total_grad.head(n_avatar_vertices() * 3).setZero();
        // writer.add_field("grad", utils::unflatten(total_grad, 3));

        // Eigen::VectorXd body_ids = Eigen::VectorXd::Zero(V.rows());
        // body_ids.head(n_avatar_vertices()).array() = 1;
        // writer.add_field("body_ids", body_ids);

        const Eigen::MatrixXd current_vertices = utils::unflatten(complete_disp, V.cols()) + V;
        // logger().debug("Save VTU to {}", path + "/step_" + std::to_string(index) + ".vtu");
        // writer.write_mesh(path + "/step_" + std::to_string(index) + ".vtu", current_vertices, F);

        garment.v = current_vertices.bottomRows(garment.v.rows());
        garment.write(path + "/step_garment_" + std::to_string(index) + ".obj");

        igl::write_triangle_mesh(path + "/step_avatar_" + std::to_string(index) + ".obj", current_vertices.topRows(n_avatar_vertices()), avatar_f);
    }

    Eigen::Vector3d bbox_size(const Eigen::Matrix<double, -1, 3> &V)
    {
        return V.colwise().maxCoeff() - V.colwise().minCoeff();
    }

    void GarmentSolver::load_garment_mesh(
        const std::string &path,
        const std::string &garment_skinning_weights_path, 
        int n_refs)
	{
        garment.read(path);

        if (std::filesystem::exists(garment_skinning_weights_path))
        {
            io::read_matrix(garment_skinning_weights_path, garment_skinning_weights);
            assert(garment_skinning_weights.rows() == skeleton_v.rows());
            assert(garment.v.rows() == garment_skinning_weights.cols());
            assert(garment_skinning_weights.minCoeff() >= 0. && garment_skinning_weights.maxCoeff() <= 1.);
        }

		// remove duplicate vertices in the garment
        // remove_duplicate_vertices(garment.v, garment.f, 1e-6);

		// while (n_refs-- > 0)
		// 	std::tie(garment.v, garment.f) = refine(garment.v, garment.f);
        assert(n_refs == 0);
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
        const std::string &target_skeleton_path,
        const std::string &target_avatar_skinning_weights_path)
    {
        igl::read_triangle_mesh(avatar_mesh_path, avatar_v, avatar_f);

        read_edge_mesh(source_skeleton_path, skeleton_v, skeleton_b);
        read_edge_mesh(target_skeleton_path, target_skeleton_v, target_skeleton_b);
        // assert((skeleton_b - target_skeleton_b).squaredNorm() < 1);
        for (int i = 0; i < skeleton_b.rows(); i++)
        {
            if ((std::min(skeleton_b(i, 0), skeleton_b(i, 1)) != std::min(target_skeleton_b(i, 0), target_skeleton_b(i, 1)))
            || (std::max(skeleton_b(i, 0), skeleton_b(i, 1)) != std::max(target_skeleton_b(i, 0), target_skeleton_b(i, 1))))
                log_and_throw_error("Inconsistent skeletons!");
        }
        target_skeleton_b = skeleton_b;

        io::read_matrix(target_avatar_skinning_weights_path, target_avatar_skinning_weights);
        assert(target_avatar_skinning_weights.rows() == skeleton_v.rows());
        assert(avatar_v.rows() == target_avatar_skinning_weights.cols());
        assert(target_avatar_skinning_weights.minCoeff() >= 0. && target_avatar_skinning_weights.maxCoeff() <= 1.);
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
        Eigen::MatrixXd new_skinning_weights(target_avatar_skinning_weights.rows(), new_avatar_v.rows());
        for (int f = 0; f < avatar_f.rows(); f++)
        {
            for (int i = 0; i < avatar_f.cols(); i++)
            {
                new_avatar_v.row(f * avatar_f.cols() + i) = avatar_v.row(avatar_f(f, i));
                new_skinning_weights.col(f * avatar_f.cols() + i) = target_avatar_skinning_weights.col(avatar_f(f, i));
                avatar_f(f, i) = f * avatar_f.cols() + i;
            }
        }
        std::swap(new_avatar_v, avatar_v);

        Eigen::VectorXi eid;
        Eigen::VectorXd coord;
        // first pass
        {
            const int N = avatar_v.rows();
            Eigen::VectorXd dists(N);
            dists.setConstant(std::numeric_limits<double>::max());
            coord.setZero(N);
            eid = -Eigen::VectorXi::Ones(N);
            for (int i = 0; i < N; i++)
            {
                Eigen::Index maxRow, maxCol;
                const double max_skin_weight = new_skinning_weights.col(i).maxCoeff(&maxRow, &maxCol);

                assert(maxCol == 0);
                for (int e = 0; e < target_skeleton_b.rows(); e++)
                {
                    if (target_skeleton_b(e, 0) == maxRow || target_skeleton_b(e, 1) == maxRow)
                    {
                        const Eigen::Vector2d tmp = project_to_edge(target_skeleton_v.row(target_skeleton_b(e, 0)), target_skeleton_v.row(target_skeleton_b(e, 1)), avatar_v.row(i));
                        if (tmp(0) < dists(i))
                        {
                            dists(i) = tmp(0);
                            coord(i) = tmp(1);
                            eid(i) = e;
                        }
                    }
                }
                if (eid(i) < 0)
                    log_and_throw_error("Failed to project vertex to the bone!");
            }

            skinny_avatar_v.setZero(avatar_v.rows(), avatar_v.cols());
            for (int i = 0; i < avatar_v.rows(); i++)
                skinny_avatar_v.row(i) += coord(i) * (skeleton_v(skeleton_b(eid(i), 1), Eigen::all) - skeleton_v(skeleton_b(eid(i), 0), Eigen::all)) + skeleton_v(skeleton_b(eid(i), 0), Eigen::all);
            skinny_avatar_f = avatar_f;
        }

        // igl::write_triangle_mesh(out_folder + "/avatar_old.obj", avatar_v, avatar_f);
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

                    if (b >= eid.size() || b < 0)
                        std::cout << "debug";

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
                        std::array<int, 2> tmp{{shared_vtx(cur, parent(source, cur)), cur}};
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
            Eigen::MatrixXd tmp_v(avatar_v.rows(), avatar_v.cols() + skinny_avatar_v.cols());
            tmp_v << avatar_v, skinny_avatar_v;
            const auto [svi, svj] = remove_duplicate_vertices(tmp_v, avatar_f, 1e-10);

            avatar_v = tmp_v.template leftCols<3>();
            skinny_avatar_v = tmp_v.template rightCols<3>();
            skinny_avatar_f = avatar_f;
        }
        
        skinny_avatar_v += (avatar_v - skinny_avatar_v) * 1e-5;
    }

    void GarmentSolver::normalize_meshes()
    {
        // Source side
        const double source_scaling = 1e2;
        skeleton_v *= source_scaling;
        garment.v *= source_scaling;
        // skinny_avatar_v *= source_scaling;

        // Target side
        const double target_scaling = bbox_size(skeleton_v).maxCoeff() / bbox_size(target_skeleton_v).maxCoeff();
        const Eigen::Vector3d center = skeleton_v.colwise().sum() / skeleton_v.rows() - target_scaling * avatar_v.colwise().sum() / avatar_v.rows();
        Transformation<3> trans(target_scaling * Eigen::Matrix3d::Identity(), center);

        trans.apply(avatar_v);
        trans.apply(target_skeleton_v);

        Eigen::MatrixXi E;
        igl::edges(garment.f, E);

        Eigen::VectorXd len = (garment.v(E.col(0), Eigen::all) - garment.v(E.col(1), Eigen::all)).rowwise().norm();
        logger().debug("Min edge length: {}", len.minCoeff());
    }
}