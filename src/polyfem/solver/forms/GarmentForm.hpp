#pragma once

#include "Form.hpp"

#include <polyfem/Common.hpp>
#include <polyfem/utils/Types.hpp>
#include <polyfem/utils/MatrixUtils.hpp>

namespace polyfem::solver
{
	class PenaltyForm : public Form
	{
	public:
		PenaltyForm(const Eigen::VectorXd &target, const std::vector<int> &indices) : target_(target), indices_(indices) { assert(indices_.size() == target_.size()); }
		virtual ~PenaltyForm() = default;

		std::string name() const override { return "Penalty"; }

	protected:
		/// @brief Compute the potential value
		/// @param x Current solution
		/// @return Value of the contact barrier potential
		double value_unweighted(const Eigen::VectorXd &x) const override
        {
            return (x(indices_) - target_).squaredNorm() / 2.;
        }

		/// @brief Compute the first derivative of the value wrt x
		/// @param[in] x Current solution
		/// @param[out] gradv Output gradient of the value wrt x
		void first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const override
        {
            gradv.setZero(x.size());
            gradv(indices_) = x(indices_) - target_;
        }

		/// @brief Compute the second derivative of the value wrt x
		/// @param x Current solution
		/// @param hessian Output Hessian of the value wrt x
		void second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const override
        {
            hessian.setZero();
			hessian.resize(x.size(), x.size());
            std::vector<Eigen::Triplet<double>> triplets;
            for (int i = 0; i < indices_.size(); i++)
                triplets.emplace_back(indices_[i], indices_[i], 1.);
            hessian.setFromTriplets(triplets.begin(), triplets.end());
        }

	protected:
        Eigen::VectorXd target_;
		std::vector<int> indices_;
	};
} // namespace polyfem::solver
