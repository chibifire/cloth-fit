#include "CurveCenterTargetForm.hpp"

#include <iostream>
#include <polyfem/utils/Logger.hpp>
#include <polyfem/autogen/auto_derivatives.hpp>
#include <polyfem/utils/AutodiffTypes.hpp>

#include <polyfem/mesh/MeshUtils.hpp>
#include <igl/write_triangle_mesh.h>
#include <finitediff.hpp>

namespace polyfem::solver
{
    namespace {
        Eigen::Vector2d point_edge_closest_distance(
            const Eigen::Vector3d &p,
            const Eigen::Vector3d &a,
            const Eigen::Vector3d &b)
        {
            const Eigen::Vector3d e = b - a;
            const Eigen::Vector3d d = p - a;
            double t = e.dot(d) / e.squaredNorm();
            t = std::min(1., std::max(0., t));
            const double dist = (d - t * e).squaredNorm();
            return Eigen::Vector2d(dist, t);
        }
        template <class T>
        Eigen::Vector<T, 2> point_line_closest_distance(
            const Eigen::Vector<T, 3> &p,
            const Eigen::Vector<T, 3> &a,
            const Eigen::Vector<T, 3> &b)
        {
            const Eigen::Vector<T, 3> e = b - a;
            const Eigen::Vector<T, 3> d = p - a;
            T t = e.dot(d) / e.squaredNorm();
            const T dist = (d - t * e).squaredNorm();
            return Eigen::Vector<T, 2>(dist, t);
        }
    }

    double OldCurveCenterTargetForm::value_unweighted(const Eigen::VectorXd &x) const 
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        double val = 0;
        for (int i = 0; i < curves_.size(); i++)
        {
            const auto &curve = curves_[i];
            const Eigen::Matrix<double, 1, 3> c = V(curve, Eigen::all).colwise().sum() / curve.size();
            val += (c - target_.row(i)).squaredNorm();
        }

        return val / 2;
    }

    void OldCurveCenterTargetForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const 
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        gradv.setZero(x.size());
        for (int i = 0; i < curves_.size(); i++)
        {
            const auto &curve = curves_[i];
            const Eigen::Matrix<double, 1, 3> c = V(curve, Eigen::all).colwise().sum() / curve.size();
            const Eigen::Vector3d g = (c - target_.row(i)) / curve.size();
            for (int j = 0; j < curve.size(); j++)
                gradv.segment(curve(j) * 3, 3) += g;
        }
    }

    void OldCurveCenterTargetForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const 
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        std::vector<Eigen::Triplet<double>> triplets;
        for (int i = 0; i < curves_.size(); i++)
        {
            const auto &curve = curves_[i];
            const double val = 1. / curve.size() / curve.size();
            for (int j = 0; j < curve.size(); j++)
                for (int d = 0; d < 3; d++)
                    for (int k = 0; k < curve.size(); k++)
                        triplets.emplace_back(curve(j) * 3 + d, curve(k) * 3 + d, val);
        }

        hessian.setZero();
        hessian.resize(x.size(), x.size());
        hessian.setFromTriplets(triplets.begin(), triplets.end());
    }

    CurveCenterProjectedTargetForm::CurveCenterProjectedTargetForm(
        const Eigen::MatrixXd &V, 
        const std::vector<Eigen::VectorXi> &curves,
        const Eigen::MatrixXd &source_skeleton_v,
        const Eigen::MatrixXd &target_skeleton_v,
        const Eigen::MatrixXi &skeleton_edges): 
        V_(V), source_skeleton_v_(source_skeleton_v),
        target_skeleton_v_(target_skeleton_v), skeleton_edges_(skeleton_edges)
    {
        for (auto curve : curves)
        {
            curves_.push_back(curve.head(curve.size()-1));
        }

        bones.resize(curves_.size());
        relative_positions.resize(curves_.size());
        for (int j = 0; j < curves_.size(); j++)
        {
            const Eigen::Vector3d center = V(curves_[j], Eigen::all).colwise().sum() / curves_[j].size();

            // Project centers to original skeleton bones
            int id = 0;
            double closest_dist = std::numeric_limits<double>::max(), closest_uv = 0;
            for (int i = 0; i < skeleton_edges.rows(); i++)
            {
                Eigen::Vector2d tmp = point_line_closest_distance<double>(center, source_skeleton_v.row(skeleton_edges(i, 0)), source_skeleton_v.row(skeleton_edges(i, 1)));
                if (tmp(0) < closest_dist)
                {
                    closest_dist = tmp(0);
                    closest_uv = tmp(1);
                    id = i;
                }
            }

            bones(j) = id;
            relative_positions(j) = closest_uv;
        }
    }

    double CurveCenterProjectedTargetForm::value_unweighted(const Eigen::VectorXd &x) const
    {
        const double t = x(0);
        const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());
        const Eigen::MatrixXd skeleton_v = source_skeleton_v_ + t * (target_skeleton_v_ - source_skeleton_v_);

        double val = 0.;
        for (int j = 0; j < curves_.size(); j++)
        {
            const Eigen::Vector3d center = V(curves_[j], Eigen::all).colwise().sum() / curves_[j].size();

            const int id = bones(j);
            const double param0 = relative_positions(j);

            const Eigen::Vector2d tmp = point_line_closest_distance<double>(center, skeleton_v.row(skeleton_edges_(id, 0)), skeleton_v.row(skeleton_edges_(id, 1)));

            val += pow(param0 - tmp(1), 2);
        }

        return val / 2.;
    }

    void CurveCenterProjectedTargetForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
    {
        const double t = x(0);
        const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());
        const Eigen::MatrixXd skeleton_v = source_skeleton_v_ + t * (target_skeleton_v_ - source_skeleton_v_);

        gradv.setZero(x.size());
        for (int j = 0; j < curves_.size(); j++)
        {
            const Eigen::Vector3d center = V(curves_[j], Eigen::all).colwise().sum() / curves_[j].size();

            const int id = bones(j);
            const double param0 = relative_positions(j);

            const Eigen::Vector2d tmp = point_line_closest_distance<double>(center, skeleton_v.row(skeleton_edges_(id, 0)), skeleton_v.row(skeleton_edges_(id, 1)));

            Eigen::Vector4d g;
            autogen::line_projection_uv_gradient(
                t, center(0), center(1), center(2), g.data(), 
                source_skeleton_v_(skeleton_edges_(id, 0), 0), source_skeleton_v_(skeleton_edges_(id, 0), 1), source_skeleton_v_(skeleton_edges_(id, 0), 2),
                target_skeleton_v_(skeleton_edges_(id, 0), 0), target_skeleton_v_(skeleton_edges_(id, 0), 1), target_skeleton_v_(skeleton_edges_(id, 0), 2),
                source_skeleton_v_(skeleton_edges_(id, 1), 0), source_skeleton_v_(skeleton_edges_(id, 1), 1), source_skeleton_v_(skeleton_edges_(id, 1), 2),
                target_skeleton_v_(skeleton_edges_(id, 1), 0), target_skeleton_v_(skeleton_edges_(id, 1), 1), target_skeleton_v_(skeleton_edges_(id, 1), 2));
            
            g *= tmp(1) - param0;

            gradv(0) += g(0);
            g.tail<3>() /= curves_[j].size();
            for (int k = 0; k < curves_[j].size(); k++)
                gradv.template segment<3>(1 + curves_[j](k) * 3) += g.tail<3>();
        }
    }

    void CurveCenterProjectedTargetForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
    {
        hessian.resize(x.size(), x.size());
        hessian.setZero();

        const double t = x(0);
        const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());
        const Eigen::MatrixXd skeleton_v = source_skeleton_v_ + t * (target_skeleton_v_ - source_skeleton_v_);

        // // debug export
        // {
        //     mesh::write_edge_mesh("skeleton.obj", skeleton_v, skeleton_edges_);
        //     Eigen::MatrixXd V1(curves_.size(), 3), V2(curves_.size(), 3);
        //     for (int j = 0; j < curves_.size(); j++)
        //     {
        //         const int N = curves_[j].size();
        //         const Eigen::Vector3d center = V(curves_[j], Eigen::all).colwise().sum() / N;

        //         const int id = bones(j);
        //         const double param0 = relative_positions(j);

        //         const Eigen::Vector2d tmp = point_line_closest_distance<double>(center, skeleton_v.row(skeleton_edges_(id, 0)), skeleton_v.row(skeleton_edges_(id, 1)));
        //         V1.row(j) = tmp(1) * (skeleton_v.row(skeleton_edges_(id, 1)) - skeleton_v.row(skeleton_edges_(id, 0))) + skeleton_v.row(skeleton_edges_(id, 0));
        //         V2.row(j) = param0 * (skeleton_v.row(skeleton_edges_(id, 1)) - skeleton_v.row(skeleton_edges_(id, 0))) + skeleton_v.row(skeleton_edges_(id, 0));
        //     }

        //     Eigen::MatrixXi F(0, 3);
        //     igl::write_triangle_mesh("current.obj", V1, F);
        //     igl::write_triangle_mesh("desired.obj", V2, F);
        // }

        std::vector<Eigen::Triplet<double>> T;
        for (int j = 0; j < curves_.size(); j++)
        {
            const auto &curve = curves_[j];
            const int N = curve.size();
            const Eigen::Vector3d center = V(curve, Eigen::all).colwise().sum() / N;

            const int id = bones(j);
            const double param0 = relative_positions(j);

            const Eigen::Vector2d tmp = point_line_closest_distance<double>(center, skeleton_v.row(skeleton_edges_(id, 0)), skeleton_v.row(skeleton_edges_(id, 1)));

            Eigen::Vector4d g;
            autogen::line_projection_uv_gradient(
                t, center(0), center(1), center(2), g.data(), 
                source_skeleton_v_(skeleton_edges_(id, 0), 0), source_skeleton_v_(skeleton_edges_(id, 0), 1), source_skeleton_v_(skeleton_edges_(id, 0), 2),
                target_skeleton_v_(skeleton_edges_(id, 0), 0), target_skeleton_v_(skeleton_edges_(id, 0), 1), target_skeleton_v_(skeleton_edges_(id, 0), 2),
                source_skeleton_v_(skeleton_edges_(id, 1), 0), source_skeleton_v_(skeleton_edges_(id, 1), 1), source_skeleton_v_(skeleton_edges_(id, 1), 2),
                target_skeleton_v_(skeleton_edges_(id, 1), 0), target_skeleton_v_(skeleton_edges_(id, 1), 1), target_skeleton_v_(skeleton_edges_(id, 1), 2));
            
            Eigen::Matrix4d h;
            autogen::line_projection_uv_hessian(
                t, center(0), center(1), center(2), h.data(), 
                source_skeleton_v_(skeleton_edges_(id, 0), 0), source_skeleton_v_(skeleton_edges_(id, 0), 1), source_skeleton_v_(skeleton_edges_(id, 0), 2),
                target_skeleton_v_(skeleton_edges_(id, 0), 0), target_skeleton_v_(skeleton_edges_(id, 0), 1), target_skeleton_v_(skeleton_edges_(id, 0), 2),
                source_skeleton_v_(skeleton_edges_(id, 1), 0), source_skeleton_v_(skeleton_edges_(id, 1), 1), source_skeleton_v_(skeleton_edges_(id, 1), 2),
                target_skeleton_v_(skeleton_edges_(id, 1), 0), target_skeleton_v_(skeleton_edges_(id, 1), 1), target_skeleton_v_(skeleton_edges_(id, 1), 2));

            h = h.eval() * (tmp(1) - param0) + g * g.transpose();

            Eigen::MatrixXd local_hess = Eigen::MatrixXd::Zero(N * 3 + 1, N * 3 + 1);
            local_hess(0, 0) = h(0, 0);
            for (int i0 = 0; i0 < N; i0++)
            {
                for (int d0 = 0; d0 < 3; d0++)
                {
                    local_hess(i0 * 3 + d0 + 1, 0) = h(d0 + 1, 0) / N;
                    local_hess(0, i0 * 3 + d0 + 1) = h(0, d0 + 1) / N;
                    for (int i1 = 0; i1 < N; i1++)
                    {
                        for (int d1 = 0; d1 < 3; d1++)
                        {
                            local_hess(i0 * 3 + d0 + 1, i1 * 3 + d1 + 1) = h(d0 + 1, d1 + 1) / (N * N);
                        }
                    }
                }
            }

            T.emplace_back(0, 0, local_hess(0, 0));
            for (int i0 = 0; i0 < N; i0++)
            for (int d0 = 0; d0 < 3; d0++)
            {
                T.emplace_back(1 + curve(i0) * 3 + d0, 0, local_hess(i0 * 3 + d0 + 1, 0));
                T.emplace_back(0, 1 + curve(i0) * 3 + d0, local_hess(0, i0 * 3 + d0 + 1));
                for (int i1 = 0; i1 < N; i1++)
                for (int d1 = 0; d1 < 3; d1++)
                {
                    T.emplace_back(1 + curve(i0) * 3 + d0, 1 + curve(i1) * 3 + d1, local_hess(i0 * 3 + d0 + 1, i1 * 3 + d1 + 1));
                }
            }
        }

        hessian.setFromTriplets(T.begin(), T.end());
    }

    CurveCenterTargetForm::CurveCenterTargetForm(
        const Eigen::MatrixXd &V, 
        const std::vector<Eigen::VectorXi> &curves,
        const Eigen::MatrixXd &source_skeleton_v,
        const Eigen::MatrixXd &target_skeleton_v,
        const Eigen::MatrixXi &skeleton_edges): 
        V_(V), source_skeleton_v_(source_skeleton_v),
        target_skeleton_v_(target_skeleton_v), skeleton_edges_(skeleton_edges)
    {
        for (auto curve : curves)
        {
            curves_.push_back(curve.head(curve.size()-1));
        }
        bones.resize(curves_.size());
        relative_positions.resize(curves_.size());
        for (int j = 0; j < curves_.size(); j++)
        {
            const Eigen::Vector3d center = V(curves_[j], Eigen::all).colwise().sum() / curves_[j].size();

            // Project centers to original skeleton bones
            int id = 0;
            double closest_dist = std::numeric_limits<double>::max(), closest_uv = 0;
            for (int i = 0; i < skeleton_edges.rows(); i++)
            {
                Eigen::Vector2d tmp = point_edge_closest_distance(center, source_skeleton_v.row(skeleton_edges(i, 0)), source_skeleton_v.row(skeleton_edges(i, 1)));
                if (tmp(0) < closest_dist)
                {
                    closest_dist = tmp(0);
                    closest_uv = tmp(1);
                    id = i;
                }
            }

            bones(j) = id;
            relative_positions(j) = closest_uv;
        }
    }

    double CurveCenterTargetForm::value_unweighted(const Eigen::VectorXd &x) const
    {
        const double t = x(0);
        const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());

        double val = 0.;
        for (int j = 0; j < curves_.size(); j++)
        {
            const Eigen::Vector3d center = V(curves_[j], Eigen::all).colwise().sum() / curves_[j].size();

            const int id = bones(j);
            const double param0 = relative_positions(j);
            const Eigen::Vector3d targetS = param0 * (source_skeleton_v_.row(skeleton_edges_(id, 1)) - source_skeleton_v_.row(skeleton_edges_(id, 0))) + source_skeleton_v_.row(skeleton_edges_(id, 0));
            const Eigen::Vector3d targetT = param0 * (target_skeleton_v_.row(skeleton_edges_(id, 1)) - target_skeleton_v_.row(skeleton_edges_(id, 0))) + target_skeleton_v_.row(skeleton_edges_(id, 0));
            const Eigen::Vector3d target = t * (targetT - targetS) + targetS;

            val += (target - center).squaredNorm();
        }

        return val / 2.;
    }

    void CurveCenterTargetForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
    {
        const double t = x(0);
        const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());

        gradv.setZero(x.size());
        for (int j = 0; j < curves_.size(); j++)
        {
            const auto &curve = curves_[j];
            const Eigen::Vector3d center = V(curve, Eigen::all).colwise().sum() / curve.size();

            const int id = bones(j);
            const double param0 = relative_positions(j);
            
            const Eigen::Vector3d targetS = param0 * (source_skeleton_v_.row(skeleton_edges_(id, 1)) - source_skeleton_v_.row(skeleton_edges_(id, 0))) + source_skeleton_v_.row(skeleton_edges_(id, 0));
            const Eigen::Vector3d targetT = param0 * (target_skeleton_v_.row(skeleton_edges_(id, 1)) - target_skeleton_v_.row(skeleton_edges_(id, 0))) + target_skeleton_v_.row(skeleton_edges_(id, 0));
            const Eigen::Vector3d target = t * (targetT - targetS) + targetS;
            
            const Eigen::Vector3d tmp = center - target;
            for (int i = 0; i < curve.size(); i++)
                gradv.segment<3>(1 + curve(i) * 3) += tmp / curve.size();
            gradv(0) -= tmp.dot(targetT - targetS);
        }
    }

    void CurveCenterTargetForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
    {
        hessian.resize(x.size(), x.size());
        hessian.setZero();

        const double t = x(0);
        const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());

        std::vector<Eigen::Triplet<double>> T;
        for (int j = 0; j < curves_.size(); j++)
        {
            const auto &curve = curves_[j];
            const int N = curve.size();
            const Eigen::Vector3d center = V(curve, Eigen::all).colwise().sum() / N;

            const int id = bones(j);
            const double param0 = relative_positions(j);
            
            const Eigen::Vector3d targetS = param0 * (source_skeleton_v_.row(skeleton_edges_(id, 1)) - source_skeleton_v_.row(skeleton_edges_(id, 0))) + source_skeleton_v_.row(skeleton_edges_(id, 0));
            const Eigen::Vector3d targetT = param0 * (target_skeleton_v_.row(skeleton_edges_(id, 1)) - target_skeleton_v_.row(skeleton_edges_(id, 0))) + target_skeleton_v_.row(skeleton_edges_(id, 0));
            const Eigen::Vector3d target = t * (targetT - targetS) + targetS;
            
            Eigen::Matrix4d h = Eigen::Matrix4d::Zero();
            h(1, 1) = 1; h(2, 2) = 1; h(3, 3) = 1;
            h(0, 0) = (targetT - targetS).dot(targetT - targetS);
            h.block<3, 1>(1, 0) = targetS - targetT;
            h.block<1, 3>(0, 1) = targetS - targetT;

            Eigen::MatrixXd local_hess = Eigen::MatrixXd::Zero(N * 3 + 1, N * 3 + 1);
            local_hess(0, 0) = h(0, 0);
            for (int i0 = 0; i0 < N; i0++)
            {
                for (int d0 = 0; d0 < 3; d0++)
                {
                    local_hess(i0 * 3 + d0 + 1, 0) = h(d0 + 1, 0) / N;
                    local_hess(0, i0 * 3 + d0 + 1) = h(0, d0 + 1) / N;
                    for (int i1 = 0; i1 < N; i1++)
                    {
                        for (int d1 = 0; d1 < 3; d1++)
                        {
                            local_hess(i0 * 3 + d0 + 1, i1 * 3 + d1 + 1) = h(d0 + 1, d1 + 1) / (N * N);
                        }
                    }
                }
            }

            // {
			// 	Eigen::MatrixXd fhess;
            //     Eigen::VectorXd x0(1 + N * 3);
            //     x0(0) = t;
            //     for (int i = 0; i < N; i++)
            //         x0.segment<3>(1 + i * 3) = V.row(curve(i));
			// 	fd::finite_hessian(
			// 		x0, [&](const Eigen::VectorXd &y) -> double 
            //         { 
            //             const Eigen::Vector3d tmp1 = y(0) * (targetT - targetS) + targetS;
            //             const Eigen::Vector3d tmp2 = utils::unflatten(y.tail(3 * N), 3).colwise().sum() / N;
            //             return (tmp1 - tmp2).squaredNorm() / 2.;
            //         }, fhess,
			// 		fd::AccuracyOrder::SECOND, 1e-6);
            //     std::cout << (fhess - local_hess).norm() / h.norm() << std::endl;
            // }
            T.emplace_back(0, 0, local_hess(0, 0));
            for (int i0 = 0; i0 < N; i0++)
            for (int d0 = 0; d0 < 3; d0++)
            {
                T.emplace_back(1 + curve(i0) * 3 + d0, 0, local_hess(i0 * 3 + d0 + 1, 0));
                T.emplace_back(0, 1 + curve(i0) * 3 + d0, local_hess(0, i0 * 3 + d0 + 1));
                for (int i1 = 0; i1 < N; i1++)
                for (int d1 = 0; d1 < 3; d1++)
                {
                    T.emplace_back(1 + curve(i0) * 3 + d0, 1 + curve(i1) * 3 + d1, local_hess(i0 * 3 + d0 + 1, i1 * 3 + d1 + 1));
                }
            }
        }

        hessian.setFromTriplets(T.begin(), T.end());
    }


    CurveTargetForm::CurveTargetForm(
        const Eigen::MatrixXd &V, 
        const std::vector<Eigen::VectorXi> &curves,
        const Eigen::MatrixXd &source_skeleton_v,
        const Eigen::MatrixXd &target_skeleton_v,
        const Eigen::MatrixXi &skeleton_edges): 
        V_(V), source_skeleton_v_(source_skeleton_v),
        target_skeleton_v_(target_skeleton_v), skeleton_edges_(skeleton_edges)
    {
        for (auto curve : curves)
        {
            curves_.push_back(curve.head(curve.size()-1));
        }
        bones.resize(curves_.size());
        relative_positions.resize(curves_.size());
        for (int j = 0; j < curves_.size(); j++)
        {
            const Eigen::Vector3d center = V(curves_[j], Eigen::all).colwise().sum() / curves_[j].size();

            // Project centers to original skeleton bones
            double closest_dist = std::numeric_limits<double>::max();
            for (int i = 0; i < skeleton_edges.rows(); i++)
            {
                Eigen::Vector2d tmp = point_edge_closest_distance(center, source_skeleton_v.row(skeleton_edges(i, 0)), source_skeleton_v.row(skeleton_edges(i, 1)));
                if (tmp(0) < closest_dist)
                {
                    closest_dist = tmp(0);
                    bones(j) = i;
                }
            }
            
            relative_positions[j].resize(curves_[j].size());
            for (int i = 0; i < curves_[j].size(); i++)
            {
                Eigen::Vector2d tmp = point_line_closest_distance<double>(V.row(curves_[j](i)), source_skeleton_v.row(skeleton_edges(bones(j), 0)), source_skeleton_v.row(skeleton_edges(bones(j), 1)));
                relative_positions[j](i) = tmp(1);
            }
        }
    }

    double CurveTargetForm::value_unweighted(const Eigen::VectorXd &x) const
    {
        const double t = x(0);
        const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());

        double val = 0.;
        for (int j = 0; j < curves_.size(); j++)
        {
            const int id = bones(j);
            for (int i = 0; i < curves_[j].size(); i++)
            {
                const Eigen::Vector3d bone0 = t * (target_skeleton_v_.row(skeleton_edges_(id, 0)) - source_skeleton_v_.row(skeleton_edges_(id, 0))) + source_skeleton_v_.row(skeleton_edges_(id, 0));
                const Eigen::Vector3d bone1 = t * (target_skeleton_v_.row(skeleton_edges_(id, 1)) - source_skeleton_v_.row(skeleton_edges_(id, 1))) + source_skeleton_v_.row(skeleton_edges_(id, 1));
                Eigen::Vector2d tmp = point_line_closest_distance<double>(V.row(curves_[j](i)), bone0, bone1);

                val += pow(relative_positions[j](i) - tmp(1), 2);
            }
        }

        return val / 2.;
    }

    void CurveTargetForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
    {
        const double t = x(0);
        const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());

		typedef DScalar1<double, Eigen::Matrix<double, 4, 1>> ADGrad;
        DiffScalarBase::setVariableCount(4);

        gradv.setZero(x.size());
        for (int j = 0; j < curves_.size(); j++)
        {
            const auto &curve = curves_[j];

            const int id = bones(j);
            for (int i = 0; i < curve.size(); i++)
            {
                Eigen::Vector<ADGrad, 4> diff_v;
                for (int d = 0; d < 3; d++)
                    diff_v(d) = ADGrad(d, V(curve(i), d));
                diff_v(3) = ADGrad(3, t);

                const Eigen::Vector3d a0 = source_skeleton_v_.row(skeleton_edges_(id, 0));
                const Eigen::Vector3d b0 = target_skeleton_v_.row(skeleton_edges_(id, 0));
                const Eigen::Vector3d a1 = source_skeleton_v_.row(skeleton_edges_(id, 1));
                const Eigen::Vector3d b1 = target_skeleton_v_.row(skeleton_edges_(id, 1));

                Eigen::Vector<ADGrad, 3> bone0;
                bone0 << diff_v(3) * (b0(0) - a0(0)) + a0(0), diff_v(3) * (b0(1) - a0(1)) + a0(1), diff_v(3) * (b0(2) - a0(2)) + a0(2);
                Eigen::Vector<ADGrad, 3> bone1;
                bone1 << diff_v(3) * (b1(0) - a1(0)) + a1(0), diff_v(3) * (b1(1) - a1(1)) + a1(1), diff_v(3) * (b1(2) - a1(2)) + a1(2);
                Eigen::Vector<ADGrad, 2> tmp = point_line_closest_distance<ADGrad>(diff_v.head<3>(), bone0, bone1);

                ADGrad err = pow(relative_positions[j](i) - tmp(1), 2) / 2.;

                gradv.segment<3>(1 + curve(i) * 3) += err.getGradient().head<3>().transpose();
                gradv(0) += err.getGradient()(3);
            }
        }
    }

    void CurveTargetForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
    {
        hessian.resize(x.size(), x.size());
        hessian.setZero();

        typedef DScalar2<double, Eigen::Matrix<double, 4, 1>, Eigen::Matrix<double, 4, 4>> ADHess;
        DiffScalarBase::setVariableCount(4);

        const double t = x(0);
        const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());

        std::vector<Eigen::Triplet<double>> T;
        for (int j = 0; j < curves_.size(); j++)
        {
            const auto &curve = curves_[j];

            const int id = bones(j);
            for (int i = 0; i < curve.size(); i++)
            {
                Eigen::Vector<ADHess, 4> diff_v;
                for (int d = 0; d < 3; d++)
                    diff_v(d) = ADHess(d, V(curve(i), d));
                diff_v(3) = ADHess(3, t);

                const Eigen::Vector3d a0 = source_skeleton_v_.row(skeleton_edges_(id, 0));
                const Eigen::Vector3d b0 = target_skeleton_v_.row(skeleton_edges_(id, 0));
                const Eigen::Vector3d a1 = source_skeleton_v_.row(skeleton_edges_(id, 1));
                const Eigen::Vector3d b1 = target_skeleton_v_.row(skeleton_edges_(id, 1));

                Eigen::Vector<ADHess, 3> bone0;
                bone0 << diff_v(3) * (b0(0) - a0(0)) + a0(0), diff_v(3) * (b0(1) - a0(1)) + a0(1), diff_v(3) * (b0(2) - a0(2)) + a0(2);
                Eigen::Vector<ADHess, 3> bone1;
                bone1 << diff_v(3) * (b1(0) - a1(0)) + a1(0), diff_v(3) * (b1(1) - a1(1)) + a1(1), diff_v(3) * (b1(2) - a1(2)) + a1(2);
                Eigen::Vector<ADHess, 2> tmp = point_line_closest_distance<ADHess>(diff_v.head<3>(), bone0, bone1);

                ADHess err = pow(relative_positions[j](i) - tmp(1), 2) / 2.;

                std::array<int, 4> indices{{curve(i) * 3 + 1, curve(i) * 3 + 2, curve(i) * 3 + 3, 0}};
                for (int k = 0; k < 4; k++)
                    for (int l = 0; l < 4; l++)
                        T.emplace_back(indices[k], indices[l], err.getHessian()(k, l));
            }
        }

        hessian.setFromTriplets(T.begin(), T.end());
    }
}