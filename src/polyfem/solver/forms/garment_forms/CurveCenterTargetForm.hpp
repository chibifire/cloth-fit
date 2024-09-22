#pragma once

#include <polyfem/solver/forms/Form.hpp>

#include <polyfem/Common.hpp>
#include <polyfem/utils/Types.hpp>
#include <polyfem/utils/MatrixUtils.hpp>

#include <vector>

namespace polyfem::solver
{
	class CurveCenterTargetForm : public Form
	{
	public:
		CurveCenterTargetForm(const Eigen::MatrixXd &V, const std::vector<Eigen::VectorXi> &curves, const Eigen::MatrixXd &target) : V_(V), curves_(curves), target_(target) {}
		virtual ~CurveCenterTargetForm() = default;

		std::string name() const override { return "curve-target"; }

	protected:
		/// @brief Compute the potential value
		/// @param x Current solution
		/// @return Value of the contact barrier potential
		double value_unweighted(const Eigen::VectorXd &x) const override;

		/// @brief Compute the first derivative of the value wrt x
		/// @param[in] x Current solution
		/// @param[out] gradv Output gradient of the value wrt x
		void first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const override;

		/// @brief Compute the second derivative of the value wrt x
		/// @param x Current solution
		/// @param hessian Output Hessian of the value wrt x
		void second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const override;

    private:
        const Eigen::MatrixXd V_;

        const std::vector<Eigen::VectorXi> curves_;
        const Eigen::MatrixXd target_;
	};
}
