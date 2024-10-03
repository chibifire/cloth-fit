#include "FitForm.hpp"

#include <openvdb/openvdb.h>
#include <openvdb/math/Vec3.h>
#include <openvdb/tools/SignedFloodFill.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/Interpolation.h>

#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/utils/Logger.hpp>

using namespace openvdb;

namespace {

    Eigen::MatrixXd upsample_standard(int N)
    {
        const int num = ((N+1)*(N+2))/2;

        Eigen::MatrixXd out(num, 4);
        for (int i = 0, k = 0; i <= N; i++)
            for (int j = 0; i + j <= N; j++)
            {
                std::array<int, 3> arr = {i, j, N - i - j};
                std::sort(arr.begin(), arr.end());

                double w;
                if (arr[1] == 0)        // vertex node
                    w = 1;
                else if (arr[0] == 0)   // edge node
                    w = 3;
                else                    // face node
                    w = 6;
                
                out.row(k) << i, j, N - i - j, w;
            }
        
        out.leftCols<3>() /= N;
        out.col(3) /= N * N;
        return out;
    }
}

namespace polyfem::solver
{
    FitForm::FitForm(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &surface_v, const Eigen::MatrixXi &surface_f, const int n_refs, const double voxel_size) : V_(V), F_(F), n_refs_(n_refs), voxel_size_(voxel_size)
    {
        math::Transform::Ptr xform = math::Transform::createLinearTransform(voxel_size);

        std::vector<Vec3s> points;
        std::vector<Vec3I> triangles;
        std::vector<Vec4I> quads;

        points.reserve(surface_v.rows());
        for (int i = 0; i < surface_v.rows(); i++)
            points.push_back(Vec3s(surface_v(i, 0), surface_v(i, 1), surface_v(i, 2)));
        
        for (int i = 0; i < surface_f.rows(); i++)
            triangles.push_back(Vec3I(surface_f(i, 0), surface_f(i, 1), surface_f(i, 2)));

        grid = tools::meshToSignedDistanceField<DoubleGrid>(*xform, points, triangles, quads, 150, 1);
        
        // build upsampling scheme on the garment mesh
        {
            Eigen::MatrixXd tmp = upsample_standard(n_refs);
            P = tmp.leftCols<3>();
            weights = tmp.col(3);
        }

        // export SDF as a triangle mesh
        // {
        //     tools::volumeToMesh(*grid, points, quads, 0.);
        //     Eigen::MatrixXd outpoints(points.size(), 3);
        //     Eigen::MatrixXi outtriangles(2 * quads.size(), 3);
        //     for (int i = 0; i < points.size(); i++)
        //     {
        //         outpoints.row(i) << points[i](0), points[i](1), points[i](2);
        //     }
        //     for (int i = 0; i < quads.size(); i++)
        //     {
        //         outtriangles.row(2 * i) << quads[i](0), quads[i](1), quads[i](2);
        //         outtriangles.row(2 * i + 1) << quads[i](0), quads[i](2), quads[i](3);
        //     }

		// 	io::OBJWriter::write(
		// 		"sdf.obj", outpoints, Eigen::MatrixXi(), outtriangles);
        // }
    }

    double FitForm::value_unweighted(const Eigen::VectorXd &x) const 
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        double val = 0;
        double max_dist = 0;
        typename DoubleGrid::ConstAccessor acc = grid->getConstAccessor();
        for (int f = 0; f < F_.rows(); f++) {
            Eigen::Matrix3d M = V({F_(f, 0),F_(f, 1),F_(f, 2)}, Eigen::all);
            Eigen::MatrixXd samples = P * M;
            const double area = (V_.row(F_(f, 1)) - V_.row(F_(f, 0))).head<3>().cross((V_.row(F_(f, 2)) - V_.row(F_(f, 0))).head<3>()).norm() / 2;
            for (int i = 0; i < P.rows(); i++) {
                math::Vec3<double> p(samples(i, 0), samples(i, 1), samples(i, 2));
                const double tmp = use_spline ? tools::SplineSampler::sample(acc, grid->transformPtr()->worldToIndex(p)) :
                                                tools::BoxSampler::sample(acc, grid->transformPtr()->worldToIndex(p));

                if (std::isnan(tmp))
                    log_and_throw_error("Invalid sdf values!");
                
                if (tmp > 0)
                    val += area * weights(i) * tmp*tmp / 2.;
            }
        }

        return val;
    }

    void FitForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const 
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        gradv.setZero(x.size());
        typename DoubleGrid::ConstAccessor acc = grid->getConstAccessor();
        for (int f = 0; f < F_.rows(); f++) {
            Eigen::Matrix3d M = V({F_(f, 0),F_(f, 1),F_(f, 2)}, Eigen::all);
            Eigen::MatrixXd samples = P * M;
            const double area = (V_.row(F_(f, 1)) - V_.row(F_(f, 0))).head<3>().cross((V_.row(F_(f, 2)) - V_.row(F_(f, 0))).head<3>()).norm() / 2;
            for (int i = 0; i < P.rows(); i++) {
                math::Vec3<double> p(samples(i, 0), samples(i, 1), samples(i, 2));
                auto tmp = use_spline ? tools::SplineSampler::sampleGradient(acc, grid->transformPtr()->worldToIndex(p)) :
                                            tools::BoxSampler::sampleGradient(acc, grid->transformPtr()->worldToIndex(p));
                tmp.g = tmp.g * (tmp.x / voxel_size_);
                
                if (tmp.x > 0)
                    for (int d = 0; d < 3; d++)
                    {
                        const double val = area * weights(i) * tmp.g(d);
                        for (int k = 0; k < P.cols(); k++)
                            gradv(3 * F_(f, k) + d) += val * P(i, k);
                    }
            }
        }
    }

    void FitForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const 
    {
        hessian.setZero();
        hessian.resize(x.size(), x.size());   
        if (!use_spline)
            return;
        
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        std::vector<Eigen::Triplet<double>> triplets;
        typename DoubleGrid::ConstAccessor acc = grid->getConstAccessor();
        for (int f = 0; f < F_.rows(); f++) {
            Eigen::Matrix3d M = V({F_(f, 0),F_(f, 1),F_(f, 2)}, Eigen::all);
            Eigen::MatrixXd samples = P * M;
            const double area = (V_.row(F_(f, 1)) - V_.row(F_(f, 0))).head<3>().cross((V_.row(F_(f, 2)) - V_.row(F_(f, 0))).head<3>()).norm() / 2;
            for (int i = 0; i < P.rows(); i++) {
                math::Vec3<double> p(samples(i, 0), samples(i, 1), samples(i, 2));
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
                h *= area * weights(i);
                if (tmp.x > 0)
                {
                    for (int d = 0; d < 3; d++)
                        for (int k = 0; k < 3; k++)
                            for (int s = 0; s < 3; s++)
                                for (int l = 0; l < 3; l++)
                                    triplets.emplace_back(F_(f, s) * 3 + d, F_(f, l) * 3 + k, P(i, s) * P(i, l) * h(d, k));
                }
            }
        }

        hessian.setFromTriplets(triplets.begin(), triplets.end());
    }
}
