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

            Eigen::VectorXi I_prime(I2.size() + 1);
            I_prime << I2, I2(1);

            curves.push_back(std::move(I_prime));
        }
        return curves;
	}

    CurveCurvatureForm::CurveCurvatureForm(const Eigen::MatrixXd &V, const std::vector<Eigen::VectorXi> &curves) : V_(V), curves_(curves)
    {
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
}
