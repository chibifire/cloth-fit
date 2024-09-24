#include "FitForm.hpp"

#include <openvdb/openvdb.h>
#include <openvdb/math/Vec3.h>
#include <openvdb/tools/SignedFloodFill.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/Interpolation.h>

#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/utils/Logger.hpp>

// #include <polyfem/utils/MaybeParallelFor.hpp>

using namespace openvdb;

// namespace
// {
//     class LocalThreadMatStorage
//     {
//     public:
//         std::unique_ptr<MatrixCache> cache = nullptr;
//         ElementAssemblyValues vals;
//         QuadratureVector da;

//         LocalThreadMatStorage() = delete;

//         LocalThreadMatStorage(const int buffer_size, const int rows, const int cols)
//         {
//             init(buffer_size, rows, cols);
//         }

//         LocalThreadMatStorage(const int buffer_size, const MatrixCache &c)
//         {
//             init(buffer_size, c);
//         }

//         LocalThreadMatStorage(const LocalThreadMatStorage &other)
//             : cache(other.cache->copy()), vals(other.vals), da(other.da)
//         {
//         }

//         LocalThreadMatStorage &operator=(const LocalThreadMatStorage &other)
//         {
//             assert(other.cache != nullptr);
//             cache = other.cache->copy();
//             vals = other.vals;
//             da = other.da;
//             return *this;
//         }

//         void init(const int buffer_size, const int rows, const int cols)
//         {
//             // assert(rows == cols);
//             // cache = std::make_unique<DenseMatrixCache>();
//             cache = std::make_unique<SparseMatrixCache>();
//             cache->reserve(buffer_size);
//             cache->init(rows, cols);
//         }

//         void init(const int buffer_size, const MatrixCache &c)
//         {
//             if (cache == nullptr)
//                 cache = c.copy();
//             cache->reserve(buffer_size);
//             cache->init(c);
//         }
//     };

//     class LocalThreadVecStorage
//     {
//     public:
//         Eigen::MatrixXd vec;
//         ElementAssemblyValues vals;
//         QuadratureVector da;

//         LocalThreadVecStorage(const int size)
//         {
//             vec.resize(size, 1);
//             vec.setZero();
//         }
//     };

//     class LocalThreadScalarStorage
//     {
//     public:
//         double val;
//         ElementAssemblyValues vals;
//         QuadratureVector da;

//         LocalThreadScalarStorage()
//         {
//             val = 0;
//         }
//     };
// } // namespace

namespace polyfem::solver
{
    FitForm::FitForm(const Eigen::MatrixXd &V, const Eigen::MatrixXd &surface_v, const Eigen::MatrixXi &surface_f, const int n_refs, const double voxel_size) : V_(V), n_refs_(n_refs), voxel_size_(voxel_size)
    {
        math::Transform::Ptr xform = math::Transform::createLinearTransform(voxel_size);

        std::vector<Vec3s> points;
        std::vector<Vec3I> triangles;
        std::vector<Vec4I> quads;

        points.reserve(surface_v.rows());
        for (int i = 0; i < surface_v.rows(); i++)
            points.push_back(Vec3s(surface_v(i, 0), surface_v(i, 1), surface_v(i, 2)) / voxel_size);
        
        for (int i = 0; i < surface_f.rows(); i++)
            triangles.push_back(Vec3I(surface_f(i, 0), surface_f(i, 1), surface_f(i, 2)));

        grid = tools::meshToSignedDistanceField<DoubleGrid>(*xform, points, triangles, quads, 50, 1);
        
        {
            tools::volumeToMesh(*grid, points, quads, 0.);
            Eigen::MatrixXd outpoints(points.size(), 3);
            Eigen::MatrixXi outtriangles(2 * quads.size(), 3);
            for (int i = 0; i < points.size(); i++)
            {
                outpoints.row(i) << points[i](0), points[i](1), points[i](2);
            }
            for (int i = 0; i < quads.size(); i++)
            {
                outtriangles.row(2 * i) << quads[i](0), quads[i](1), quads[i](2);
                outtriangles.row(2 * i + 1) << quads[i](0), quads[i](2), quads[i](3);
            }

			io::OBJWriter::write(
				"sdf.obj", outpoints, Eigen::MatrixXi(), outtriangles);
        }
    }

    double FitForm::value_unweighted(const Eigen::VectorXd &x) const 
    {
        Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        double val = 0;
        double max_dist = 0;
        for (int i = 0; i < V.rows(); i++) {
            typename DoubleGrid::ConstAccessor acc = grid->getConstAccessor();
            math::Vec3<double> p(V(i, 0), V(i, 1), V(i, 2));
            const double tmp = use_spline ? tools::SplineSampler::sample(acc, grid->transformPtr()->worldToIndex(p)) :
                                            tools::BoxSampler::sample(acc, grid->transformPtr()->worldToIndex(p));

            if (std::isnan(tmp))
                log_and_throw_error("Invalid sdf values!");
            
            if (tmp > 0)
                val += tmp*tmp / 2.;
        }

        return val;
    }

    void FitForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const 
    {
        Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        gradv.setZero(x.size());
        for (int i = 0; i < V.rows(); i++) {
            typename DoubleGrid::ConstAccessor acc = grid->getConstAccessor();
            math::Vec3<double> p(V(i, 0), V(i, 1), V(i, 2));
            auto tmp = use_spline ? tools::SplineSampler::sampleGradient(acc, grid->transformPtr()->worldToIndex(p)) :
                                          tools::BoxSampler::sampleGradient(acc, grid->transformPtr()->worldToIndex(p));
            tmp.g = tmp.g * (tmp.x / voxel_size_);
            
            if (tmp.x > 0)
                for (int d = 0; d < 3; d++)
                    gradv(3 * i + d) += tmp.g(d);
        }
    }

    void FitForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const 
    {
        hessian.setZero();
        hessian.resize(x.size(), x.size());   
        if (!use_spline)
            return;
        
        Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        std::vector<Eigen::Triplet<double>> triplets;
        for (int i = 0; i < V.rows(); i++) {
            typename DoubleGrid::ConstAccessor acc = grid->getConstAccessor();
            math::Vec3<double> p(V(i, 0), V(i, 1), V(i, 2));
            auto tmp = tools::SplineSampler::sampleHessian(acc, grid->transformPtr()->worldToIndex(p));
            tmp.g = tmp.g / voxel_size_;
            
            Eigen::Vector3d g;
            g << tmp.g(0), tmp.g(1), tmp.g(2);
            Eigen::Matrix3d h;
            h << tmp.h(0, 0), tmp.h(0, 1), tmp.h(0, 2), 
                tmp.h(1, 0), tmp.h(1, 1), tmp.h(1, 2),
                tmp.h(2, 0), tmp.h(2, 1), tmp.h(2, 2);
            h *= tmp.x / voxel_size_ / voxel_size_;
            h += g * g.transpose();
            if (tmp.x > 0)
                for (int d = 0; d < 3; d++)
                    for (int k = 0; k < 3; k++)
                        triplets.emplace_back(i * 3 + d, i * 3 + k, h(d, k));
        }

        hessian.setFromTriplets(triplets.begin(), triplets.end());
    }
}
