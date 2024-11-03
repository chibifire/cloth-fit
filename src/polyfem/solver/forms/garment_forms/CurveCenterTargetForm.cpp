#include "CurveCenterTargetForm.hpp"

#include <iostream>
#include <polyfem/utils/Logger.hpp>
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
        Eigen::Vector2d point_line_closest_distance(
            const Eigen::Vector3d &p,
            const Eigen::Vector3d &a,
            const Eigen::Vector3d &b)
        {
            const Eigen::Vector3d e = b - a;
            const Eigen::Vector3d d = p - a;
            double t = e.dot(d) / e.squaredNorm();
            const double dist = (d - t * e).squaredNorm();
            return Eigen::Vector2d(dist, t);
        }
        void line_projection_uv_gradient(
            double t, double vx, double vy, double vz, double dA[4],
            double a0x, double a0y, double a0z,
            double a1x, double a1y, double a1z,
            double b0x, double b0y, double b0z,
            double b1x, double b1y, double b1z)
        {
            const auto t0 = -b0x;
            const auto t1 = b0x - b1x;
            const auto t2 = a0x - a1x;
            const auto t3 = a0x - t * t2;
            const auto t4 = t * t1 + t0 + t3;
            const auto t5 = -b0y;
            const auto t6 = b0y - b1y;
            const auto t7 = a0y - a1y;
            const auto t8 = a0y - t * t7;
            const auto t9 = t * t6 + t5 + t8;
            const auto t10 = -b0z;
            const auto t11 = b0z - b1z;
            const auto t12 = a0z - a1z;
            const auto t13 = a0z - t * t12;
            const auto t14 = t * t11 + t10 + t13;
            const auto t15 =
                1.0 / (std::pow(t14, 2) + std::pow(t4, 2) + std::pow(t9, 2));
            const auto t16 = b1x + t0 + t2;
            const auto t17 = t3 - vx;
            const auto t18 = b1y + t5 + t7;
            const auto t19 = t8 - vy;
            const auto t20 = b1z + t10 + t12;
            const auto t21 = t13 - vz;
            const auto t22 = a0x + t * t1 - t * t2 + t0;
            const auto t23 = a0y + t * t6 - t * t7 + t5;
            const auto t24 = a0z + t * t11 - t * t12 + t10;
            const auto t25 =
                1.0 / (std::pow(t22, 2) + std::pow(t23, 2) + std::pow(t24, 2));
            dA[0] = -t15
                * (t12 * t14
                + 2 * t15 * (-t14 * t20 - t16 * t4 - t18 * t9)
                    * (t14 * t21 + t17 * t4 + t19 * t9)
                + t16 * t17 + t18 * t19 + t2 * t4 + t20 * t21 + t7 * t9);
            dA[1] = -t22 * t25;
            dA[2] = -t23 * t25;
            dA[3] = -t24 * t25;
        }

        // dA is (16×1) flattened in column-major order
        void line_projection_uv_hessian(
            double t, double vx, double vy, double vz, double dA[16],
            double a0x, double a0y, double a0z,
            double a1x, double a1y, double a1z,
            double b0x, double b0y, double b0z,
            double b1x, double b1y, double b1z)
        {
            const auto t0 = a0x - a1x;
            const auto t1 = -b0x;
            const auto t2 = b1x + t0 + t1;
            const auto t3 = a0y - a1y;
            const auto t4 = -b0y;
            const auto t5 = b1y + t3 + t4;
            const auto t6 = a0z - a1z;
            const auto t7 = -b0z;
            const auto t8 = b1z + t6 + t7;
            const auto t9 = b0x - b1x;
            const auto t10 = t * t0;
            const auto t11 = a0x + t * t9 + t1 - t10;
            const auto t12 = b0y - b1y;
            const auto t13 = t * t3;
            const auto t14 = a0y + t * t12 - t13 + t4;
            const auto t15 = b0z - b1z;
            const auto t16 = t * t6;
            const auto t17 = a0z + t * t15 - t16 + t7;
            const auto t18 = std::pow(t11, 2) + std::pow(t14, 2) + std::pow(t17, 2);
            const auto t19 = 1.0 / t18;
            const auto t20 = -a0x;
            const auto t21 = -a0y;
            const auto t22 = -a0z;
            const auto t23 = t11 * (-t10 - t20 - vx) + t14 * (-t13 - t21 - vy)
                + t17 * (-t16 - t22 - vz);
            const auto t24 = a0x - t * t0;
            const auto t25 = t * t9 + t1 + t24;
            const auto t26 = a0y - t * t3;
            const auto t27 = t * t12 + t26 + t4;
            const auto t28 = a0z - t * t6;
            const auto t29 = t * t15 + t28 + t7;
            const auto t30 = t2 * t25 + t27 * t5 + t29 * t8;
            const auto t31 = 2 * t19;
            const auto t32 = t30 * t31;
            const auto t33 = t19 * (-a1x - t11 * t32 - t20 - t9);
            const auto t34 = t19 * (-a1y - t12 - t14 * t32 - t21);
            const auto t35 = t19 * (-a1z - t15 - t17 * t32 - t22);
            dA[0] = t31
                * (t0 * t2
                - t19 * t23
                    * (std::pow(t2, 2) + std::pow(t5, 2) + std::pow(t8, 2))
                + t3 * t5
                - t32
                    * (t0 * t25 + t2 * (t24 - vx) + t27 * t3 + t29 * t6
                        + t5 * (t26 - vy) + t8 * (t28 - vz))
                + t6 * t8 + 4 * t23 * std::pow(t30, 2) / std::pow(t18, 2));
            dA[1] = t33;
            dA[2] = t34;
            dA[3] = t35;
            dA[4] = t33;
            dA[5] = 0;
            dA[6] = 0;
            dA[7] = 0;
            dA[8] = t34;
            dA[9] = 0;
            dA[10] = 0;
            dA[11] = 0;
            dA[12] = t35;
            dA[13] = 0;
            dA[14] = 0;
            dA[15] = 0;
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
                Eigen::Vector2d tmp = point_line_closest_distance(center, source_skeleton_v.row(skeleton_edges(i, 0)), source_skeleton_v.row(skeleton_edges(i, 1)));
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

            const Eigen::Vector2d tmp = point_line_closest_distance(center, skeleton_v.row(skeleton_edges_(id, 0)), skeleton_v.row(skeleton_edges_(id, 1)));

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

            const Eigen::Vector2d tmp = point_line_closest_distance(center, skeleton_v.row(skeleton_edges_(id, 0)), skeleton_v.row(skeleton_edges_(id, 1)));

            Eigen::Vector4d g;
            line_projection_uv_gradient(
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

        //         const Eigen::Vector2d tmp = point_line_closest_distance(center, skeleton_v.row(skeleton_edges_(id, 0)), skeleton_v.row(skeleton_edges_(id, 1)));
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

            const Eigen::Vector2d tmp = point_line_closest_distance(center, skeleton_v.row(skeleton_edges_(id, 0)), skeleton_v.row(skeleton_edges_(id, 1)));

            Eigen::Vector4d g;
            line_projection_uv_gradient(
                t, center(0), center(1), center(2), g.data(), 
                source_skeleton_v_(skeleton_edges_(id, 0), 0), source_skeleton_v_(skeleton_edges_(id, 0), 1), source_skeleton_v_(skeleton_edges_(id, 0), 2),
                target_skeleton_v_(skeleton_edges_(id, 0), 0), target_skeleton_v_(skeleton_edges_(id, 0), 1), target_skeleton_v_(skeleton_edges_(id, 0), 2),
                source_skeleton_v_(skeleton_edges_(id, 1), 0), source_skeleton_v_(skeleton_edges_(id, 1), 1), source_skeleton_v_(skeleton_edges_(id, 1), 2),
                target_skeleton_v_(skeleton_edges_(id, 1), 0), target_skeleton_v_(skeleton_edges_(id, 1), 1), target_skeleton_v_(skeleton_edges_(id, 1), 2));
            
            Eigen::Matrix4d h;
            line_projection_uv_hessian(
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
                gradv.segment<3>(1 + curve(i) * 3, 3) += tmp / curve.size();
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
}