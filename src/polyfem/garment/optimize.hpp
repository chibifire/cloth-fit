#pragma once
#include <polyfem/Common.hpp>

namespace ipc {
    class CollisionMesh;
}

namespace polyfem {
    namespace solver {
        class GarmentNLProblem;
    }

    void save_vtu(
        const std::string &path,
        solver::GarmentNLProblem &prob,
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const int n_avatar_vertices,
        const Eigen::VectorXd &sol);
    
    Eigen::Vector3d bbox_size(const Eigen::Matrix<double, -1, 3> &V);

    class GarmentSolver 
    {
    public:
        void check_intersections(
            const ipc::CollisionMesh &collision_mesh,
            const Eigen::MatrixXd &collision_vertices) const;

        void load_garment_mesh(
            const std::string &path,
            int n_refs = 0);

        void read_meshes(
            const std::string &avatar_mesh_path,
            const std::string &source_skeleton_path,
            const std::string &target_skeleton_path);
        
        void project_avatar_to_skeleton();

        void normalize_meshes();
    
        std::string out_folder;
    
        Eigen::MatrixXd avatar_v, garment_v;
        Eigen::MatrixXi avatar_f, garment_f;

        Eigen::MatrixXd skeleton_v, target_skeleton_v;
        Eigen::MatrixXi skeleton_b, target_skeleton_b;

        Eigen::MatrixXd skinny_avatar_v;
        Eigen::MatrixXi skinny_avatar_f;
    };
}