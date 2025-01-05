#pragma once
#include <polyfem/Common.hpp>

namespace ipc {
    class CollisionMesh;
}

namespace polyfem {
    namespace solver {
        class GarmentNLProblem;
    }

    struct OBJMesh {
        // vertices, triangle faces
        Eigen::MatrixXd v;
        Eigen::MatrixXi f;
        // normals
        Eigen::MatrixXi fn;
        Eigen::MatrixXd cn;
        // textures
        Eigen::MatrixXi ftc;
        Eigen::MatrixXd tc;

        void read(const std::string &path);
        void write(const std::string &path);
    };
    
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
            const std::string &target_skeleton_path,
            const std::string &target_avatar_skinning_weights_path);
        
        void project_avatar_to_skeleton();

        void normalize_meshes();
    
        void save_result(
            const std::string &path,
            const int index,
            solver::GarmentNLProblem &prob,
            const Eigen::MatrixXd &V,
            const Eigen::MatrixXi &F,
            const Eigen::VectorXd &sol);
        
        int n_garment_vertices() const { return garment.v.rows(); }
        int n_garment_faces() const { return garment.f.rows(); }
        int n_avatar_vertices() const { return avatar_v.rows(); }

        std::string out_folder;
    
        Eigen::MatrixXd avatar_v;
        Eigen::MatrixXi avatar_f;

        OBJMesh garment;

        Eigen::MatrixXd skeleton_v, target_skeleton_v;
        Eigen::MatrixXi skeleton_b, target_skeleton_b;

        Eigen::MatrixXd skinny_avatar_v;
        Eigen::MatrixXi skinny_avatar_f;

        Eigen::MatrixXd target_avatar_skinning_weights;
        Eigen::MatrixXd garment_skinning_weights;

        std::vector<int> not_fit_fids;
    };
}