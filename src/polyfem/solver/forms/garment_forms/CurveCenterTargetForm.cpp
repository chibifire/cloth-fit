#include "CurveCenterTargetForm.hpp"

#include <polyfem/utils/Logger.hpp>

namespace polyfem::solver
{
    double CurveCenterTargetForm::value_unweighted(const Eigen::VectorXd &x) const 
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

    void CurveCenterTargetForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const 
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

    void CurveCenterTargetForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const 
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

}