#include "CurveConstraintForm.hpp"

#include <polyfem/utils/Logger.hpp>
#include <polyfem/autogen/auto_derivatives.hpp>

#include <igl/boundary_facets.h>
#include <igl/adjacency_matrix.h>
#include <igl/connected_components.h>
#include <igl/edges_to_path.h>

using namespace polyfem::autogen;

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
    }
	std::vector<Eigen::VectorXi> boundary_curves(const Eigen::MatrixXi &F)
	{
        Eigen::MatrixXi edges;
        igl::boundary_facets(F, edges);

        if (edges.size() == 0)
            return {};
        
        Eigen::SparseMatrix<int> A;
        igl::adjacency_matrix(edges, A);

        Eigen::VectorXi C, K1;
        const int n_connected = igl::connected_components(A, C, K1);

        // logger().debug("[{}] Number of curve loops: {}", name(), n_connected);

        std::vector<Eigen::VectorXi> curves;
        for (int i = 0; i < n_connected; i++)
        {
            if (K1(i) < 2)
                continue;

            std::vector<int> ind;            
            for (int j = 0; j < edges.rows(); j++)
            {
                if (C(edges(j, 0)) == i)
                    ind.push_back(j);
            }

            const Eigen::MatrixXi edges_tmp = edges(ind, Eigen::all);

            Eigen::VectorXi I2, J2, K2;
            igl::edges_to_path(edges_tmp, I2, J2, K2);
            assert(I2(0) == I2(I2.size() - 1));
            assert(I2.size() >= 4);

            // Eigen::VectorXi I_prime(I2.size() + 1);
            // I_prime << I2, I2(1);

            curves.push_back(std::move(I2));
        }
        return curves;
	}

    Eigen::MatrixXd extract_curve_center_targets(
        const Eigen::MatrixXd &garment_v,
        const std::vector<Eigen::VectorXi> &curves,
        const Eigen::MatrixXd &skeleton_v,
        const Eigen::MatrixXi &skeleton_bones,
        const Eigen::MatrixXd &target_skeleton_v)
    {
        Eigen::MatrixXd targets(curves.size(), 3);
        for (int j = 0; j < curves.size(); j++)
        {
            // Compute centers of curves
            Eigen::Vector3d center = garment_v(curves[j], Eigen::all).colwise().sum() / curves[j].size();

            // Project centers to original skeleton bones
            int id = 0;
            double closest_dist = std::numeric_limits<double>::max(), closest_uv = 0;
            for (int i = 0; i < skeleton_bones.rows(); i++)
            {
                Eigen::Vector2d tmp = point_edge_closest_distance(center, skeleton_v.row(skeleton_bones(i, 0)), skeleton_v.row(skeleton_bones(i, 1)));
                if (tmp(0) < closest_dist)
                {
                    closest_dist = tmp(0);
                    closest_uv = tmp(1);
                    id = i;
                }
            }

            // Map positions to new skeleton bones
            targets.row(j) = closest_uv * (target_skeleton_v.row(skeleton_bones(id, 1)) - target_skeleton_v.row(skeleton_bones(id, 0))) + target_skeleton_v.row(skeleton_bones(id, 0));
        }

        return targets;
    }

    CurveCurvatureForm::CurveCurvatureForm(const Eigen::MatrixXd &V, const std::vector<Eigen::VectorXi> &curves) : V_(V)
    {
        for (const auto &c : curves)
        {
            assert(c(0) == c(c.size() - 1));

            Eigen::VectorXi c_(c.size() + 1);
            c_.head(c.size()) = c;
            c_(c.size()) = c(1);

            assert(c_(0) == c_(c_.size() - 2));
            assert(c_(1) == c_(c_.size() - 1));

            curves_.push_back(c_);
        }
        orig_angles = compute_angles(V);
    }

    CurveCurvatureForm::CurveCurvatureForm(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) : CurveCurvatureForm(V, boundary_curves(F))
    {
    }

    std::vector<Eigen::MatrixXd> CurveCurvatureForm::compute_angles(const Eigen::MatrixXd &V) const
    {
        std::vector<Eigen::MatrixXd> out;
        for (const auto &curve : curves_)
        {
            Eigen::MatrixXd result(curve.size() - 2, 2);
            for (int i = 0; i < curve.size() - 2; i++)
            {
                const Eigen::Ref<const Eigen::Vector3d> a = V.row(curve(i));
                const Eigen::Ref<const Eigen::Vector3d> b = V.row(curve(i+1));
                const Eigen::Ref<const Eigen::Vector3d> c = V.row(curve(i+2));
                const Eigen::Vector3d v1 = c - b;
                const Eigen::Vector3d v2 = a - b;
                const double len = v1.squaredNorm() * v2.squaredNorm();

                result(i, 0) = pow(v1.dot(v2), 2) / len;
                result(i, 1) = v1.cross(v2).squaredNorm() / len;
            }
            out.push_back(std::move(result));
        }

        return out;
    }

    double CurveCurvatureForm::value_unweighted(const Eigen::VectorXd &x) const
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;
        auto angles = compute_angles(V);

        double val = 0;
        for (int i = 0; i < angles.size(); i++)
            val += (angles[i] - orig_angles[i]).squaredNorm() / 2.;
        
        return val;
    }

    void CurveCurvatureForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        gradv.setZero(x.size());
        for (int c = 0; c < curves_.size(); c++)
        {
            const auto &curve = curves_[c];
            for (int i = 0; i < curve.size() - 2; i++)
            {
                const Eigen::Ref<const Eigen::Vector3d> p1 = V.row(curve(i));
                const Eigen::Ref<const Eigen::Vector3d> p2 = V.row(curve(i+1));
                const Eigen::Ref<const Eigen::Vector3d> p3 = V.row(curve(i+2));
                const Eigen::Vector3d v1 = p3 - p2;
                const Eigen::Vector3d v2 = p1 - p2;
                const double len = v1.squaredNorm() * v2.squaredNorm();
                const double errA = pow(v1.dot(v2), 2) / len - orig_angles[c](i, 0);
                const double errB = v1.cross(v2).squaredNorm() / len - orig_angles[c](i, 1);

                Eigen::Vector<double, 9> g;
                curve_dot_product_norm_gradient(
                    p1(0), p1(1), p1(2),
                    p2(0), p2(1), p2(2),
                    p3(0), p3(1), p3(2), g.data());
                
                for (int k = 0; k < 3; k++)
                    gradv.segment(curve(i + k) * 3, 3) += g.segment<3>(k * 3) * errA;

                curve_cross_product_norm_gradient(
                    p1(0), p1(1), p1(2),
                    p2(0), p2(1), p2(2),
                    p3(0), p3(1), p3(2), g.data());

                for (int k = 0; k < 3; k++)
                    gradv.segment(curve(i + k) * 3, 3) += g.segment<3>(k * 3) * errB;
            }
        }
    }

    void CurveCurvatureForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        std::vector<Eigen::Triplet<double>> triplets;
        for (int c = 0; c < curves_.size(); c++)
        {
            const auto &curve = curves_[c];
            for (int i = 0; i < curve.size() - 2; i++)
            {
                const Eigen::Ref<const Eigen::Vector3d> p1 = V.row(curve(i));
                const Eigen::Ref<const Eigen::Vector3d> p2 = V.row(curve(i+1));
                const Eigen::Ref<const Eigen::Vector3d> p3 = V.row(curve(i+2));
                const Eigen::Vector3d v1 = p3 - p2;
                const Eigen::Vector3d v2 = p1 - p2;
                const double len = v1.squaredNorm() * v2.squaredNorm();
                const double errA = pow(v1.dot(v2), 2) / len - orig_angles[c](i, 0);
                const double errB = v1.cross(v2).squaredNorm() / len - orig_angles[c](i, 1);

                Eigen::Vector<double, 9> g;
                Eigen::Matrix<double, 9, 9> h, local_hess;

                curve_dot_product_norm_gradient(
                    p1(0), p1(1), p1(2),
                    p2(0), p2(1), p2(2),
                    p3(0), p3(1), p3(2), g.data());

                curve_dot_product_norm_hessian(
                    p1(0), p1(1), p1(2),
                    p2(0), p2(1), p2(2),
                    p3(0), p3(1), p3(2), h.data());
                
                local_hess = h * errA + g * g.transpose();

                curve_cross_product_norm_gradient(
                    p1(0), p1(1), p1(2),
                    p2(0), p2(1), p2(2),
                    p3(0), p3(1), p3(2), g.data());

                curve_cross_product_norm_hessian(
                    p1(0), p1(1), p1(2),
                    p2(0), p2(1), p2(2),
                    p3(0), p3(1), p3(2), h.data());
                
                local_hess += h * errB + g * g.transpose();

                for (int lj = 0; lj < 3; lj++)
                    for (int dj = 0; dj < 3; dj++)
                        for (int li = 0; li < 3; li++)
                            for (int di = 0; di < 3; di++)
                                triplets.emplace_back(curve(i + li) * 3 + di, curve(i + lj) * 3 + dj, local_hess(li * 3 + di, lj * 3 + dj));
            }
        }
    
        hessian.setZero();
        hessian.resize(x.size(), x.size());
        hessian.setFromTriplets(triplets.begin(), triplets.end());
    }


    CurveTwistForm::CurveTwistForm(const Eigen::MatrixXd &V, const std::vector<Eigen::VectorXi> &curves) : V_(V)
    {
        for (const auto &c : curves)
        {
            const int N = c.size() - 1;
            assert(c(0) == c(N));

            Eigen::VectorXi c_(N + 3);
            c_.head(N + 1) = c;
            c_(N + 1) = c(1);
            c_(N + 2) = c(2);

            assert(c_(0) == c_(c_.size() - 3));
            assert(c_(1) == c_(c_.size() - 2));
            assert(c_(2) == c_(c_.size() - 1));

            curves_.push_back(c_);
        }
        orig_angles = compute_angles(V);
    }

    CurveTwistForm::CurveTwistForm(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) : CurveTwistForm(V, boundary_curves(F))
    {
    }

    std::vector<Eigen::VectorXd> CurveTwistForm::compute_angles(const Eigen::MatrixXd &V) const
    {
        std::vector<Eigen::VectorXd> out;
        for (const auto &curve : curves_)
        {
            Eigen::VectorXd result(curve.size() - 3);
            for (int i = 0; i < curve.size() - 3; i++)
            {
                const Eigen::Ref<const Eigen::Vector3d> a = V.row(curve(i));
                const Eigen::Ref<const Eigen::Vector3d> b = V.row(curve(i+1));
                const Eigen::Ref<const Eigen::Vector3d> c = V.row(curve(i+2));
                const Eigen::Ref<const Eigen::Vector3d> d = V.row(curve(i+3));
                const Eigen::Vector3d mid = (c - b).normalized();
                const Eigen::Vector3d v1 = b - a;
                const Eigen::Vector3d v2 = d - c;

                const Eigen::Vector3d w1 = (v1 - v1.dot(mid) * mid).normalized();
                const Eigen::Vector3d w2 = (v2 - v2.dot(mid) * mid).normalized();

                result(i) = 1. - w1.dot(w2);
            }
            out.push_back(std::move(result));
        }

        return out;
    }

    double CurveTwistForm::value_unweighted(const Eigen::VectorXd &x) const
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;
        auto angles = compute_angles(V);

        double val = 0;
        for (int i = 0; i < angles.size(); i++)
            val += (angles[i] - orig_angles[i]).squaredNorm() / 2.;
        
        return val;
    }

    void CurveTwistForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        gradv.setZero(x.size());
        for (int c = 0; c < curves_.size(); c++)
        {
            const auto &curve = curves_[c];
            for (int i = 0; i < curve.size() - 3; i++)
            {
                const Eigen::Ref<const Eigen::Vector3d> p0 = V.row(curve(i));
                const Eigen::Ref<const Eigen::Vector3d> p1 = V.row(curve(i+1));
                const Eigen::Ref<const Eigen::Vector3d> p2 = V.row(curve(i+2));
                const Eigen::Ref<const Eigen::Vector3d> p3 = V.row(curve(i+3));
                const Eigen::Vector3d mid = (p2 - p1).normalized();
                const Eigen::Vector3d v1 = p1 - p0;
                const Eigen::Vector3d v2 = p3 - p2;
                const Eigen::Vector3d w1 = (v1 - v1.dot(mid) * mid).normalized();
                const Eigen::Vector3d w2 = (v2 - v2.dot(mid) * mid).normalized();

                const double err = (1. - w1.dot(w2)) - orig_angles[c](i);

                Eigen::Vector<double, 12> g;
                curve_twist_gradient(
                    p0(0), p0(1), p0(2),
                    p1(0), p1(1), p1(2),
                    p2(0), p2(1), p2(2), 
                    p3(0), p3(1), p3(2), g.data());
                
                if (!std::isfinite(g.squaredNorm()))
                    log_and_throw_error("NAN curve twist");

                for (int k = 0; k < 4; k++)
                    gradv.segment(curve(i + k) * 3, 3) += g.segment<3>(k * 3) * err;
            }
        }
    }

    void CurveTwistForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        std::vector<Eigen::Triplet<double>> triplets;
        for (int c = 0; c < curves_.size(); c++)
        {
            const auto &curve = curves_[c];
            for (int i = 0; i < curve.size() - 3; i++)
            {
                const Eigen::Ref<const Eigen::Vector3d> p0 = V.row(curve(i));
                const Eigen::Ref<const Eigen::Vector3d> p1 = V.row(curve(i+1));
                const Eigen::Ref<const Eigen::Vector3d> p2 = V.row(curve(i+2));
                const Eigen::Ref<const Eigen::Vector3d> p3 = V.row(curve(i+3));
                const Eigen::Vector3d mid = (p2 - p1).normalized();
                const Eigen::Vector3d v1 = p1 - p0;
                const Eigen::Vector3d v2 = p3 - p2;
                const Eigen::Vector3d w1 = (v1 - v1.dot(mid) * mid).normalized();
                const Eigen::Vector3d w2 = (v2 - v2.dot(mid) * mid).normalized();

                const double err = (1. - w1.dot(w2)) - orig_angles[c](i);

                Eigen::Vector<double, 12> g;
                Eigen::Matrix<double, 12, 12> h;

                curve_twist_gradient(
                    p0(0), p0(1), p0(2),
                    p1(0), p1(1), p1(2),
                    p2(0), p2(1), p2(2), 
                    p3(0), p3(1), p3(2), g.data());

                curve_twist_hessian(
                    p0(0), p0(1), p0(2),
                    p1(0), p1(1), p1(2),
                    p2(0), p2(1), p2(2), 
                    p3(0), p3(1), p3(2), h.data());
                
                h = (h * err + g * g.transpose()).eval();

                for (int lj = 0; lj < 4; lj++)
                    for (int dj = 0; dj < 3; dj++)
                        for (int li = 0; li < 4; li++)
                            for (int di = 0; di < 3; di++)
                                triplets.emplace_back(curve(i + li) * 3 + di, curve(i + lj) * 3 + dj, h(li * 3 + di, lj * 3 + dj));
            }
        }
    
        hessian.setZero();
        hessian.resize(x.size(), x.size());
        hessian.setFromTriplets(triplets.begin(), triplets.end());
    }

    SymmetryForm::SymmetryForm(const Eigen::MatrixXd &V, const Eigen::VectorXi &curve): V_(V)
    {
        {
            assert(curve(0) == curve(curve.size() - 1));
            curve_ = curve.head(curve.size() - 1);
        }

        const Eigen::MatrixXd P = V(curve_, Eigen::all);
        const Eigen::Vector3d bbox_min = P.colwise().minCoeff();
        const Eigen::Vector3d bbox_max = P.colwise().maxCoeff();
        const double bbox_size = (bbox_max - bbox_min).maxCoeff();

        const Eigen::Vector3d center = V(curve_, Eigen::all).colwise().sum() / curve_.size();

        std::vector<Eigen::Vector3d> coordinates;
        for (int i = 0; i < curve_.size(); i++)
        {
            Eigen::Vector3d x = P.row(i);
            x(dim) = P(i, dim) - center(dim);
            coordinates.push_back(x);
        }
        
        correspondence.setZero(coordinates.size());
        double max_err = 0;
        for (int i = 0; i < coordinates.size(); i++)
        {
            Eigen::Vector3d x = coordinates[i];
            bool found = false;
            double min_err = std::numeric_limits<double>::max();
            int min_id = -1;
            for (int j = 0; j < coordinates.size(); j++)
            {
                Eigen::Vector3d y = coordinates[j];
                double err = abs(x(dim) + y(dim));
                for (int d = 0; d < 3; d++)
                    if (d != dim)
                        err += abs(x(d) - y(d));
                if (err < tol * bbox_size)
                {
                    max_err = std::max(max_err, err);
                    found = true;
                    correspondence(i) = j;
                    break;
                }
                if (err < min_err)
                {
                    min_err = err;
                    min_id = j;
                }
            }
            if (!found)
            {
                logger().error("Asymmetric vertex (ID {}, pos {}) on the curve with error {} (ID {}, pos {}) found! Set weight to zero!", 
                    curve_(i), x.transpose(), min_err / bbox_size, curve_(min_id), coordinates[min_id].transpose());

                disable();
                break;
            }
        }

        if (enabled())
        {
            logger().debug("Symmetric curve identified! Error is {}", max_err / bbox_size);
        }
    }

    double SymmetryForm::value_unweighted(const Eigen::VectorXd &x) const
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;
        const Eigen::Vector3d center = V(curve_, Eigen::all).colwise().sum() / curve_.size();

        Eigen::MatrixXd tmp = V(curve_(correspondence), Eigen::all) - V(curve_, Eigen::all);
        tmp.col(dim) = (V(curve_(correspondence), dim) + V(curve_, dim)).array() - 2 * center(dim);

        return tmp.squaredNorm() / 2.;
    }

    void SymmetryForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const 
    {
        const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;
        const Eigen::Vector3d center = V(curve_, Eigen::all).colwise().sum() / curve_.size();

        Eigen::MatrixXd tmp = V(curve_(correspondence), Eigen::all) - V(curve_, Eigen::all);
        tmp.col(dim) = (V(curve_(correspondence), dim) + V(curve_, dim)).array() - 2 * center(dim);

        Eigen::MatrixXd g = Eigen::MatrixXd::Zero(V.rows(), V.cols());
        for (int d = 0; d < 3; d++)
        {
            if (d != dim)
            {
                g(curve_(correspondence), d) += tmp.col(d);
                g(curve_, d) -= tmp.col(d);
            }
        }
        
        g(curve_(correspondence), dim) += tmp.col(dim);
        g(curve_, dim) += tmp.col(dim);

        const double deriv_wrt_c = -2 * tmp.col(dim).sum();
        g(curve_, dim).array() += deriv_wrt_c / curve_.size();

        gradv = utils::flatten(g);
    }

    void SymmetryForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
    {
        const int N = curve_.size();
        // const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

        hessian.setZero();
        hessian.resize(x.size(), x.size());
        std::vector<Eigen::Triplet<double>> T;
        
        Eigen::MatrixXd H(N, N);
        for (int d = 0; d < 3; d++)
        {
            H.setZero();
            H.diagonal().array() += 2;
            
            for (int i = 0; i < N; i++)
                H(i, correspondence(i)) += 2 * (d == dim ? 1 : -1);
            
            if (d == dim)
                H.array() += (-4. / N);

            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    if (H(i, j) != 0)
                        T.emplace_back(curve_(i) * 3 + d, curve_(j) * 3 + d, H(i, j));
        }

        hessian.setFromTriplets(T.begin(), T.end());
    }
}
