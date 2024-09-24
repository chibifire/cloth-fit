#pragma once

#include <polyfem/solver/forms/Form.hpp>

#include <polyfem/Common.hpp>
#include <polyfem/utils/Types.hpp>
#include <polyfem/utils/MatrixUtils.hpp>

#include <openvdb/openvdb.h>

namespace polyfem::solver
{
	class FitForm : public Form
	{
	public:
		FitForm(const Eigen::MatrixXd &V, const Eigen::MatrixXd &surface_v, const Eigen::MatrixXi &surface_f, const int n_refs, const double voxel_size);

		std::string name() const override { return "garment-fit"; }

		/// @brief Compute the value of the form
		/// @param x Current solution
		/// @return Computed value
		double value_unweighted(const Eigen::VectorXd &x) const override;

		/// @brief Compute the first derivative of the value wrt x
		/// @param[in] x Current solution
		/// @param[out] gradv Output gradient of the value wrt x
		void first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const override;

		/// @brief Compute the second derivative of the value wrt x
		/// @param[in] x Current solution
		/// @param[out] hessian Output Hessian of the value wrt x
		void second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const override;

    private:
        const Eigen::MatrixXd V_;
        const int n_refs_;
        const double voxel_size_;
        const bool use_spline = true;

        openvdb::DoubleGrid::Ptr grid;
	};
} // namespace polyfem::solver
