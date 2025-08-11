#include "GarmentForm.hpp"
#include <Eigen/Dense>
#include <igl/triangle_triangle_adjacency.h>
#include <polyfem/utils/Logger.hpp>
#include <polyfem/utils/Timer.hpp>
#include <polyfem/autogen/auto_derivatives.hpp>
#include <unordered_set>

using namespace polyfem::autogen;

namespace
{
	Eigen::Vector3d triangle_normal(const Eigen::RowVector3d &a, const Eigen::RowVector3d &b, const Eigen::RowVector3d &c)
	{
		return (b - a).cross(c - a);
	}

	Eigen::MatrixXd compute_normals(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
	{
		Eigen::MatrixXd normals = Eigen::MatrixXd::Zero(F.rows(), 3);
		for (int f = 0; f < F.rows(); f++)
		{
			normals.row(f) = triangle_normal(V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)));
		}

		return normals;
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

	Eigen::Vector3d vector_in_affine_coordinate(
		const Eigen::Vector3d &a,
		const Eigen::Vector3d &b,
		const Eigen::Vector3d &c,
		const Eigen::Vector3d &v)
	{
		Eigen::Matrix3d A;
		A << a, b, c;
		Eigen::Vector3d x = A.colPivHouseholderQr().solve(v);
		return x;
	}
} // namespace

namespace polyfem::solver
{
	double AreaForm::value_unweighted(const Eigen::VectorXd &x) const
	{
		double total = 0;
		const Eigen::MatrixXd V_new = utils::unflatten(x, V_.cols()) + V_;
		for (int e = 0; e < F_.rows(); e++)
		{
			const Eigen::Vector3d a = V_new.row(F_(e, 1)) - V_new.row(F_(e, 0));
			const Eigen::Vector3d b = V_new.row(F_(e, 2)) - V_new.row(F_(e, 0));
			const double area = a.cross(b).norm() / 2.;
			const double x = area / threshold_;
			if (x < 1.)
				total -= std::log(x) * (1. - x) * (1. - x);
		}
		return total;
	}

	void AreaForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	{
		gradv.setZero(x.size());
		const Eigen::MatrixXd V_new = utils::unflatten(x, V_.cols()) + V_;
		for (int e = 0; e < F_.rows(); e++)
		{
			const Eigen::Vector3d a = V_new.row(F_(e, 1)) - V_new.row(F_(e, 0));
			const Eigen::Vector3d b = V_new.row(F_(e, 2)) - V_new.row(F_(e, 0));
			const double area = a.cross(b).norm() / 2.;
			const double x = area / threshold_;
			if (x < 1.)
			{
				Eigen::Vector<double, 9> local_grad;
				local_grad.setZero();
				triangle_area_gradient(
					V_new(F_(e, 0), 0), V_new(F_(e, 0), 1), V_new(F_(e, 0), 2),
					V_new(F_(e, 1), 0), V_new(F_(e, 1), 1), V_new(F_(e, 1), 2),
					V_new(F_(e, 2), 0), V_new(F_(e, 2), 1), V_new(F_(e, 2), 2),
					local_grad.data());
				// derivative of barrier wrt. area
				local_grad *= -(x - 1) * (2 * x * std::log(x) + x - 1) / area;
				for (int i = 0; i < 3; i++)
					gradv.segment(F_(e, i) * 3, 3) += local_grad.segment(i * 3, 3);
			}
		}
	}

	void AreaForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
	{
		POLYFEM_SCOPED_TIMER("area hessian");
		hessian.setZero();
		hessian.resize(x.size(), x.size());
		std::vector<Eigen::Triplet<double>> triplets;

		const Eigen::MatrixXd V_new = utils::unflatten(x, V_.cols()) + V_;
		for (int e = 0; e < F_.rows(); e++)
		{
			const Eigen::Vector3d a = V_new.row(F_(e, 1)) - V_new.row(F_(e, 0));
			const Eigen::Vector3d b = V_new.row(F_(e, 2)) - V_new.row(F_(e, 0));
			const double area = a.cross(b).norm() / 2.;
			const double x = area / threshold_;
			if (x < 1.)
			{
				Eigen::Vector<double, 9> local_grad;
				local_grad.setZero();
				triangle_area_gradient(
					V_new(F_(e, 0), 0), V_new(F_(e, 0), 1), V_new(F_(e, 0), 2),
					V_new(F_(e, 1), 0), V_new(F_(e, 1), 1), V_new(F_(e, 1), 2),
					V_new(F_(e, 2), 0), V_new(F_(e, 2), 1), V_new(F_(e, 2), 2),
					local_grad.data());

				Eigen::Matrix<double, 9, 9> local_hess;
				local_hess.setZero();
				triangle_area_hessian(
					V_new(F_(e, 0), 0), V_new(F_(e, 0), 1), V_new(F_(e, 0), 2),
					V_new(F_(e, 1), 0), V_new(F_(e, 1), 1), V_new(F_(e, 1), 2),
					V_new(F_(e, 2), 0), V_new(F_(e, 2), 1), V_new(F_(e, 2), 2),
					local_hess.data());

				// derivative of barrier wrt. area
				const double barrier_grad = -(x - 1) * (2 * x * std::log(x) + x - 1) / area;
				const double barrier_hess = -(2 * std::log(x) + 3 - 2 / x - 1 / x / x) / threshold_;

				Eigen::Matrix<double, 9, 9> real_local_hess = local_grad * (local_grad.transpose() * barrier_hess) + local_hess * barrier_grad;
				if (is_project_to_psd())
					real_local_hess = ipc::project_to_psd(real_local_hess);
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						for (int di = 0; di < 3; di++)
							for (int dj = 0; dj < 3; dj++)
								triplets.emplace_back(3 * F_(e, i) + di, 3 * F_(e, j) + dj, real_local_hess(i * 3 + di, j * 3 + dj));
			}
		}

		hessian.setFromTriplets(triplets.begin(), triplets.end());
	}

	double DefGradForm::value_unweighted(const Eigen::VectorXd &x) const
	{
		double total = 0;
		const Eigen::MatrixXd V = utils::unflatten(x, V_.cols()) + V_;
		for (int e = 0; e < F_.rows(); e++)
		{
			double a0 = (V_.row(F_(e, 1)) - V_.row(F_(e, 0))).norm();
			double b0 = (V_.row(F_(e, 2)) - V_.row(F_(e, 0))).norm();
			double a1 = (V.row(F_(e, 1)) - V.row(F_(e, 0))).norm();
			double b1 = (V.row(F_(e, 2)) - V.row(F_(e, 0))).norm();

			total += pow(a1 / b1 - a0 / b0, 2);
		}
		return total / 2.;
	}

	void DefGradForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	{
		gradv.setZero(x.size());
		const Eigen::MatrixXd V = utils::unflatten(x, V_.cols()) + V_;
		for (int e = 0; e < F_.rows(); e++)
		{
			double a0 = (V_.row(F_(e, 1)) - V_.row(F_(e, 0))).norm();
			double b0 = (V_.row(F_(e, 2)) - V_.row(F_(e, 0))).norm();
			double a1 = (V.row(F_(e, 1)) - V.row(F_(e, 0))).norm();
			double b1 = (V.row(F_(e, 2)) - V.row(F_(e, 0))).norm();

			Eigen::Matrix<double, 9, 1> local_grad;
			def_grad_gradient(
				V(F_(e, 0), 0), V(F_(e, 0), 1), V(F_(e, 0), 2),
				V(F_(e, 1), 0), V(F_(e, 1), 1), V(F_(e, 1), 2),
				V(F_(e, 2), 0), V(F_(e, 2), 1), V(F_(e, 2), 2),
				local_grad.data());

			for (int li = 0; li < 3; li++)
				gradv.segment<3>(F_(e, li) * 3) += (a1 / b1 - a0 / b0) * local_grad.segment<3>(li * 3);
		}
	}

	void DefGradForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
	{
		POLYFEM_SCOPED_TIMER("def grad hessian");
		hessian.setZero();
		hessian.resize(x.size(), x.size());
		std::vector<Eigen::Triplet<double>> triplets;

		const Eigen::MatrixXd V = utils::unflatten(x, V_.cols()) + V_;
		for (int e = 0; e < F_.rows(); e++)
		{
			double a0 = (V_.row(F_(e, 1)) - V_.row(F_(e, 0))).norm();
			double b0 = (V_.row(F_(e, 2)) - V_.row(F_(e, 0))).norm();
			double a1 = (V.row(F_(e, 1)) - V.row(F_(e, 0))).norm();
			double b1 = (V.row(F_(e, 2)) - V.row(F_(e, 0))).norm();

			Eigen::Vector<double, 9> local_grad;
			local_grad.setZero();
			def_grad_gradient(
				V(F_(e, 0), 0), V(F_(e, 0), 1), V(F_(e, 0), 2),
				V(F_(e, 1), 0), V(F_(e, 1), 1), V(F_(e, 1), 2),
				V(F_(e, 2), 0), V(F_(e, 2), 1), V(F_(e, 2), 2),
				local_grad.data());

			Eigen::Matrix<double, 9, 9> local_hess;
			local_hess.setZero();
			def_grad_hessian(
				V(F_(e, 0), 0), V(F_(e, 0), 1), V(F_(e, 0), 2),
				V(F_(e, 1), 0), V(F_(e, 1), 1), V(F_(e, 1), 2),
				V(F_(e, 2), 0), V(F_(e, 2), 1), V(F_(e, 2), 2),
				local_hess.data());

			const Eigen::Matrix<double, 9, 9> real_local_hess = local_grad * local_grad.transpose() + local_hess * (a1 / b1 - a0 / b0);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					for (int di = 0; di < 3; di++)
						for (int dj = 0; dj < 3; dj++)
							triplets.emplace_back(3 * F_(e, i) + di, 3 * F_(e, j) + dj, real_local_hess(i * 3 + di, j * 3 + dj));
		}

		hessian.setFromTriplets(triplets.begin(), triplets.end());
	}

	NormalForm::NormalForm(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) : V_(V), F_(F)
	{
		Eigen::MatrixXd normals = compute_normals(V, F_);
		orig_areas = normals.rowwise().norm();
		orig_areas /= 2;
	}

	double NormalForm::value_unweighted(const Eigen::VectorXd &x) const
	{
		double total = 0;
		const Eigen::MatrixXd V = utils::unflatten(x, V_.cols()) + V_;
		for (int e = 0; e < F_.rows(); e++)
		{
			Eigen::Vector3d a0 = V_.row(F_(e, 1)) - V_.row(F_(e, 0));
			Eigen::Vector3d b0 = V_.row(F_(e, 2)) - V_.row(F_(e, 0));
			Eigen::Vector3d n0 = a0.cross(b0).normalized();
			Eigen::Vector3d a1 = V.row(F_(e, 1)) - V.row(F_(e, 0));
			Eigen::Vector3d b1 = V.row(F_(e, 2)) - V.row(F_(e, 0));
			Eigen::Vector3d n1 = a1.cross(b1).normalized();

			total += orig_areas(e) * (1 - n0.dot(n1));
		}
		return total;
	}

	void NormalForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	{
		gradv.setZero(x.size());
		const Eigen::MatrixXd V = utils::unflatten(x, V_.cols()) + V_;
		for (int e = 0; e < F_.rows(); e++)
		{
			Eigen::Vector3d a0 = V_.row(F_(e, 1)) - V_.row(F_(e, 0));
			Eigen::Vector3d b0 = V_.row(F_(e, 2)) - V_.row(F_(e, 0));
			Eigen::Vector3d n0 = a0.cross(b0).normalized();

			Eigen::Matrix<double, 12, 1> local_grad;
			normal_gradient(
				V(F_(e, 0), 0), V(F_(e, 0), 1), V(F_(e, 0), 2),
				V(F_(e, 1), 0), V(F_(e, 1), 1), V(F_(e, 1), 2),
				V(F_(e, 2), 0), V(F_(e, 2), 1), V(F_(e, 2), 2),
				n0(0), n0(1), n0(2),
				local_grad.data());

			for (int li = 0; li < 3; li++)
				gradv.segment<3>(F_(e, li) * 3) += orig_areas(e) * local_grad.segment<3>(li * 3);
		}
	}

	void NormalForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
	{
		POLYFEM_SCOPED_TIMER("normal hessian");
		hessian.setZero();
		hessian.resize(x.size(), x.size());
		std::vector<Eigen::Triplet<double>> triplets;

		const Eigen::MatrixXd V = utils::unflatten(x, V_.cols()) + V_;
		for (int e = 0; e < F_.rows(); e++)
		{
			Eigen::Vector3d a0 = V_.row(F_(e, 1)) - V_.row(F_(e, 0));
			Eigen::Vector3d b0 = V_.row(F_(e, 2)) - V_.row(F_(e, 0));
			Eigen::Vector3d n0 = a0.cross(b0).normalized();

			Eigen::Matrix<double, 12, 12> local_hess;
			local_hess.setZero();
			normal_hessian(
				V(F_(e, 0), 0), V(F_(e, 0), 1), V(F_(e, 0), 2),
				V(F_(e, 1), 0), V(F_(e, 1), 1), V(F_(e, 1), 2),
				V(F_(e, 2), 0), V(F_(e, 2), 1), V(F_(e, 2), 2),
				n0(0), n0(1), n0(2),
				local_hess.data());

			if (is_project_to_psd())
				local_hess = ipc::project_to_psd(local_hess);

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					for (int di = 0; di < 3; di++)
						for (int dj = 0; dj < 3; dj++)
							triplets.emplace_back(3 * F_(e, i) + di, 3 * F_(e, j) + dj, orig_areas(e) * local_hess(i * 3 + di, j * 3 + dj));
		}

		hessian.setFromTriplets(triplets.begin(), triplets.end());
	}

	AngleForm::AngleForm(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) : V_(V), F_(F)
	{
		igl::triangle_triangle_adjacency(F_, TT, TTi);

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;
				// std::unordered_set<int> edge1 = {F_(TT(i, j), le(TTi(i, j), 0)), F_(TT(i, j), le(TTi(i, j), 1))};
				std::unordered_set<int> edge2 = {F_(i, le(j, 0)), F_(i, le(j, 1))};
				if (TTi(i, j) < 0)
				{
					for (int k = 0; k < 3; k++)
					{
						std::unordered_set<int> edge = {F_(TT(i, j), le(k, 0)), F_(TT(i, j), le(k, 1))};
						if (edge == edge2)
						{
							TTi(i, j) = k;
							break;
						}
					}
				}
			}
		}

		Eigen::MatrixXd normals = compute_normals(V, F_);

		areas = normals.rowwise().norm() / 2.;

		orig_angles = Eigen::MatrixXd::Ones(TT.rows(), TT.cols() * 2);
		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;

				const double normal_len = normals.row(i).norm() * normals.row(TT(i, j)).norm();
				orig_angles(i, 2 * j) = normals.row(i).dot(normals.row(TT(i, j))) / normal_len;
				const Eigen::Vector3d shared_edge = V.row(F_(i, le(j, 1))) - V.row(F_(i, le(j, 0)));
				orig_angles(i, 2 * j + 1) = normals.row(i).head<3>().cross(normals.row(TT(i, j)).head<3>()).dot(shared_edge) / shared_edge.norm() / normal_len;
			}
		}
	}

	double AngleForm::value_unweighted(const Eigen::VectorXd &x) const
	{
		Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;
		Eigen::MatrixXd normals = compute_normals(V, F_);

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		double result = 0;
		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;

				const double normal_len = normals.row(i).norm() * normals.row(TT(i, j)).norm();
				const double errA = normals.row(i).dot(normals.row(TT(i, j))) / normal_len - orig_angles(i, 2 * j);
				const Eigen::Vector3d shared_edge = V.row(F_(i, le(j, 1))) - V.row(F_(i, le(j, 0)));
				const double errB = normals.row(i).head<3>().cross(normals.row(TT(i, j)).head<3>()).dot(shared_edge) / shared_edge.norm() / normal_len - orig_angles(i, 2 * j + 1);

				result += (errA * errA + errB * errB) * (areas(i) + areas(TT(i, j)));
			}
		}
		return result / 2.;
	}

	void AngleForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	{
		Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;
		Eigen::MatrixXd normals = compute_normals(V, F_);

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		gradv.setZero(V.size());
		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;

				assert(TTi(i, j) >= 0);

				const double normal_len_i = normals.row(i).norm();
				const double normal_len_j = normals.row(TT(i, j)).norm();
				const double normal_len = normal_len_i * normal_len_j;
				const double dot_prod = normals.row(i).dot(normals.row(TT(i, j)));
				const double errA = dot_prod / normal_len - orig_angles(i, 2 * j);
				const Eigen::Vector3d shared_edge = V.row(F_(i, le(j, 1))) - V.row(F_(i, le(j, 0)));
				const double errB = normals.row(i).head<3>().cross(normals.row(TT(i, j)).head<3>()).dot(shared_edge) / shared_edge.norm() / normal_len - orig_angles(i, 2 * j + 1);
				const double fac = areas(i) + areas(TT(i, j));

				Eigen::Vector4i indices;
				assert(TTi(i, j) < lv.size());
				assert(TT(i, j) < F_.rows());
				indices << F_(i, le(j, 0)), F_(i, le(j, 1)), F_(i, lv(j)), F_(TT(i, j), lv(TTi(i, j)));
				assert(F_(TT(i, j), lv(TTi(i, j))) != F_(i, le(j, 0)));
				assert(F_(TT(i, j), lv(TTi(i, j))) != F_(i, le(j, 1)));

				// grad of cosine
				{
					Eigen::Matrix<double, 12, 1> local_grad;
					normal_dot_product_gradient(
						V(indices(0), 0), V(indices(0), 1), V(indices(0), 2),
						V(indices(1), 0), V(indices(1), 1), V(indices(1), 2),
						V(indices(2), 0), V(indices(2), 1), V(indices(2), 2),
						V(indices(3), 0), V(indices(3), 1), V(indices(3), 2),
						local_grad.data());

					for (int li = 0; li < indices.size(); li++)
						gradv.segment<3>(indices(li) * 3) += fac * errA * local_grad.segment<3>(li * 3);
				}

				// grad of sine
				{
					Eigen::Matrix<double, 12, 1> local_grad;

					normal_triple_product_gradient(
						V(indices(0), 0), V(indices(0), 1), V(indices(0), 2),
						V(indices(1), 0), V(indices(1), 1), V(indices(1), 2),
						V(indices(2), 0), V(indices(2), 1), V(indices(2), 2),
						V(indices(3), 0), V(indices(3), 1), V(indices(3), 2),
						local_grad.data());

					for (int li = 0; li < indices.size(); li++)
						gradv.segment<3>(indices(li) * 3) += fac * errB * local_grad.segment<3>(li * 3);
				}
			}
		}
	}

	void AngleForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
	{
		POLYFEM_SCOPED_TIMER("angle hessian");
		Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;
		Eigen::MatrixXd normals = compute_normals(V, F_);

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		std::vector<Eigen::Triplet<double>> triplets;
		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;

				const double normal_len_i = normals.row(i).norm();
				const double normal_len_j = normals.row(TT(i, j)).norm();
				const double normal_len = normal_len_i * normal_len_j;
				const double dot_prod = normals.row(i).dot(normals.row(TT(i, j)));
				const double errA = dot_prod / normal_len - orig_angles(i, 2 * j);
				const Eigen::Vector3d shared_edge = V.row(F_(i, le(j, 1))) - V.row(F_(i, le(j, 0)));
				const double errB = normals.row(i).head<3>().cross(normals.row(TT(i, j)).head<3>()).dot(shared_edge) / shared_edge.norm() / normal_len - orig_angles(i, 2 * j + 1);
				const double fac = areas(i) + areas(TT(i, j));

				Eigen::Vector4i indices;
				indices << F_(i, le(j, 0)), F_(i, le(j, 1)), F_(i, lv(j)), F_(TT(i, j), lv(TTi(i, j)));
				assert(F_(TT(i, j), lv(TTi(i, j))) != F_(i, le(j, 0)));
				assert(F_(TT(i, j), lv(TTi(i, j))) != F_(i, le(j, 1)));

				Eigen::Matrix<double, 12, 12> local_hess;
				local_hess.setZero();

				// hess of cosine
				{
					Eigen::Matrix<double, 12, 1> local_grad;
					normal_dot_product_gradient(
						V(indices(0), 0), V(indices(0), 1), V(indices(0), 2),
						V(indices(1), 0), V(indices(1), 1), V(indices(1), 2),
						V(indices(2), 0), V(indices(2), 1), V(indices(2), 2),
						V(indices(3), 0), V(indices(3), 1), V(indices(3), 2),
						local_grad.data());
					Eigen::Matrix<double, 12, 12> tmp_hess;
					normal_dot_product_hessian(
						V(indices(0), 0), V(indices(0), 1), V(indices(0), 2),
						V(indices(1), 0), V(indices(1), 1), V(indices(1), 2),
						V(indices(2), 0), V(indices(2), 1), V(indices(2), 2),
						V(indices(3), 0), V(indices(3), 1), V(indices(3), 2),
						tmp_hess.data());

					local_hess += fac * (tmp_hess * errA + local_grad * local_grad.transpose());
				}

				// grad of sine
				{
					Eigen::Matrix<double, 12, 1> local_grad;
					normal_triple_product_gradient(
						V(indices(0), 0), V(indices(0), 1), V(indices(0), 2),
						V(indices(1), 0), V(indices(1), 1), V(indices(1), 2),
						V(indices(2), 0), V(indices(2), 1), V(indices(2), 2),
						V(indices(3), 0), V(indices(3), 1), V(indices(3), 2),
						local_grad.data());

					Eigen::Matrix<double, 12, 12> tmp_hess;
					normal_triple_product_hessian(
						V(indices(0), 0), V(indices(0), 1), V(indices(0), 2),
						V(indices(1), 0), V(indices(1), 1), V(indices(1), 2),
						V(indices(2), 0), V(indices(2), 1), V(indices(2), 2),
						V(indices(3), 0), V(indices(3), 1), V(indices(3), 2),
						tmp_hess.data());

					local_hess += fac * (tmp_hess * errB + local_grad * local_grad.transpose());
				}

				if (is_project_to_psd())
					local_hess = ipc::project_to_psd(local_hess);

				for (int lj = 0; lj < indices.size(); lj++)
					for (int dj = 0; dj < 3; dj++)
						for (int li = 0; li < indices.size(); li++)
							for (int di = 0; di < 3; di++)
								triplets.emplace_back(indices(li) * 3 + di, indices(lj) * 3 + dj, local_hess(li * 3 + di, lj * 3 + dj));
			}
		}

		hessian.setZero();
		hessian.resize(x.size(), x.size());
		hessian.setFromTriplets(triplets.begin(), triplets.end());
	}

	RelativeScalingForm::RelativeScalingForm(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) : V_(V), F_(F)
	{
		igl::triangle_triangle_adjacency(F_, TT, TTi);

		{
			Eigen::MatrixXd normals = compute_normals(V, F_);
			orig_areas = normals.rowwise().norm();
			orig_areas /= 2;
		}

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;
				// std::unordered_set<int> edge1 = {F_(TT(i, j), le(TTi(i, j), 0)), F_(TT(i, j), le(TTi(i, j), 1))};
				std::unordered_set<int> edge2 = {F_(i, le(j, 0)), F_(i, le(j, 1))};
				if (TTi(i, j) < 0)
				{
					for (int k = 0; k < 3; k++)
					{
						std::unordered_set<int> edge = {F_(TT(i, j), le(k, 0)), F_(TT(i, j), le(k, 1))};
						if (edge == edge2)
						{
							TTi(i, j) = k;
							break;
						}
					}
				}
			}
		}

		orig_dists = Eigen::MatrixXd::Ones(TT.rows(), TT.cols() * 2);
		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;

				const double edge_len = (V.row(F_(i, le(j, 1))) - V.row(F_(i, le(j, 0)))).norm();
				orig_dists(i, j * 2 + 0) = orig_areas(i) / edge_len;
				orig_dists(i, j * 2 + 1) = orig_areas(TT(i, j)) / edge_len;
			}
		}
	}

	double RelativeScalingForm::value_unweighted(const Eigen::VectorXd &x) const
	{
		const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

		Eigen::VectorXd areas;
		{
			Eigen::MatrixXd normals = compute_normals(V, F_);
			areas = normals.rowwise().norm();
		}

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		double result = 0;
		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;

				const double edge_len = (V.row(F_(i, le(j, 1))) - V.row(F_(i, le(j, 0)))).norm();
				const double err = areas(i) / edge_len / orig_dists(i, 2 * j) - areas(TT(i, j)) / edge_len / orig_dists(i, 2 * j + 1);
				result += err * err * (orig_areas(i) + orig_areas(TT(i, j)));
			}
		}

		return result / 2.;
	}

	void RelativeScalingForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	{
		const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

		Eigen::VectorXd areas;
		{
			Eigen::MatrixXd normals = compute_normals(V, F_);
			areas = normals.rowwise().norm();
		}

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		gradv.setZero(V.size());
		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;

				const double edge_len = (V.row(F_(i, le(j, 1))) - V.row(F_(i, le(j, 0)))).norm();
				const double err = areas(i) / edge_len / orig_dists(i, 2 * j) - areas(TT(i, j)) / edge_len / orig_dists(i, 2 * j + 1);
				const double fac = orig_areas(i) + orig_areas(TT(i, j));

				Eigen::Matrix<int, 2, 3> indices;
				indices << F_(i, le(j, 0)), F_(i, le(j, 1)), F_(i, lv(j)),
					F_(i, le(j, 0)), F_(i, le(j, 1)), F_(TT(i, j), lv(TTi(i, j)));

				for (int k = 0; k < 2; k++)
				{
					Eigen::Matrix<double, 9, 1> local_grad;
					point_edge_distance_gradient(
						V(indices(k, 0), 0), V(indices(k, 0), 1), V(indices(k, 0), 2),
						V(indices(k, 1), 0), V(indices(k, 1), 1), V(indices(k, 1), 2),
						V(indices(k, 2), 0), V(indices(k, 2), 1), V(indices(k, 2), 2),
						local_grad.data());

					for (int li = 0; li < indices.cols(); li++)
						gradv.segment<3>(indices(k, li) * 3) += (fac * ((k == 0) ? 1 : -1) * err / orig_dists(i, 2 * j + k)) * local_grad.segment<3>(li * 3);
				}
			}
		}
	}

	void RelativeScalingForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
	{
		POLYFEM_SCOPED_TIMER("similarity hessian");
		const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

		Eigen::VectorXd areas;
		{
			Eigen::MatrixXd normals = compute_normals(V, F_);
			areas = normals.rowwise().norm();
		}

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		std::vector<Eigen::Triplet<double>> triplets;
		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;

				const double edge_len = (V.row(F_(i, le(j, 1))) - V.row(F_(i, le(j, 0)))).norm();
				const double err = areas(i) / edge_len / orig_dists(i, 2 * j) - areas(TT(i, j)) / edge_len / orig_dists(i, 2 * j + 1);
				const double fac = orig_areas(i) + orig_areas(TT(i, j));

				Eigen::Matrix<int, 2, 3> indices;
				indices << F_(i, le(j, 0)), F_(i, le(j, 1)), F_(i, lv(j)),
					F_(i, le(j, 0)), F_(i, le(j, 1)), F_(TT(i, j), lv(TTi(i, j)));

				std::array<Eigen::Matrix<double, 9, 1>, 2> local_grad;
				std::array<Eigen::Matrix<double, 9, 9>, 2> local_hess;
				for (int k = 0; k < 2; k++)
				{
					point_edge_distance_gradient(
						V(indices(k, 0), 0), V(indices(k, 0), 1), V(indices(k, 0), 2),
						V(indices(k, 1), 0), V(indices(k, 1), 1), V(indices(k, 1), 2),
						V(indices(k, 2), 0), V(indices(k, 2), 1), V(indices(k, 2), 2),
						local_grad[k].data());

					point_edge_distance_hessian(
						V(indices(k, 0), 0), V(indices(k, 0), 1), V(indices(k, 0), 2),
						V(indices(k, 1), 0), V(indices(k, 1), 1), V(indices(k, 1), 2),
						V(indices(k, 2), 0), V(indices(k, 2), 1), V(indices(k, 2), 2),
						local_hess[k].data());
				}

				Eigen::Vector<double, 12> g;
				g.setZero();
				g({0, 1, 2, 3, 4, 5, 6, 7, 8}) += (1. / orig_dists(i, 2 * j + 0)) * local_grad[0];
				g({0, 1, 2, 3, 4, 5, 9, 10, 11}) -= (1. / orig_dists(i, 2 * j + 1)) * local_grad[1];
				Eigen::Matrix<double, 12, 12> h = (fac * g) * g.transpose();

				h({0, 1, 2, 3, 4, 5, 6, 7, 8}, {0, 1, 2, 3, 4, 5, 6, 7, 8}) += (fac * err / orig_dists(i, 2 * j + 0)) * local_hess[0];
				h({0, 1, 2, 3, 4, 5, 9, 10, 11}, {0, 1, 2, 3, 4, 5, 9, 10, 11}) -= (fac * err / orig_dists(i, 2 * j + 1)) * local_hess[1];

				if (is_project_to_psd())
					h = ipc::project_to_psd(h);

				Eigen::Vector<int, 4> indices4;
				indices4 << F_(i, le(j, 0)), F_(i, le(j, 1)), F_(i, lv(j)), F_(TT(i, j), lv(TTi(i, j)));
				for (int li = 0; li < indices4.size(); li++)
					for (int di = 0; di < 3; di++)
						for (int lj = 0; lj < indices4.size(); lj++)
							for (int dj = 0; dj < 3; dj++)
								if (h(li * 3 + di, lj * 3 + dj) != 0)
									triplets.emplace_back(3 * indices4(li) + di, 3 * indices4(lj) + dj, h(li * 3 + di, lj * 3 + dj));
			}
		}

		hessian.setZero();
		hessian.resize(x.size(), x.size());
		hessian.setFromTriplets(triplets.begin(), triplets.end());
	}

	SimilarityForm::SimilarityForm(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) : V_(V), F_(F)
	{
		igl::triangle_triangle_adjacency(F_, TT, TTi);

		{
			Eigen::MatrixXd normals = compute_normals(V, F_);
			orig_areas = normals.rowwise().norm();
			orig_areas /= 2;
		}

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		for (int i = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++)
			{
				if (TT(i, j) < 0)
					continue;
				// std::unordered_set<int> edge1 = {F_(TT(i, j), le(TTi(i, j), 0)), F_(TT(i, j), le(TTi(i, j), 1))};
				std::unordered_set<int> edge2 = {F_(i, le(j, 0)), F_(i, le(j, 1))};
				if (TTi(i, j) < 0)
				{
					for (int k = 0; k < 3; k++)
					{
						std::unordered_set<int> edge = {F_(TT(i, j), le(k, 0)), F_(TT(i, j), le(k, 1))};
						if (edge == edge2)
						{
							TTi(i, j) = k;
							break;
						}
					}
				}
			}
		}

		orig_coeffs.setZero(TT.size(), 6);
		for (int i = 0, k = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++, k++)
			{
				if (TT(i, j) < 0)
					continue;

				const Eigen::Vector3d v = (V.row(F_(TT(i, j), lv(TTi(i, j)))) - V.row(F_(i, lv(j)))).normalized();
				const Eigen::Vector3d e = (V.row(F_(i, le(j, 1))) - V.row(F_(i, le(j, 0))));
				const Eigen::Vector3d e0 = (V.row(F_(i, lv(j))) - V.row(F_(i, le(j, 0))));
				const Eigen::Vector3d e1 = (V.row(F_(TT(i, j), lv(TTi(i, j)))) - V.row(F_(i, le(j, 0))));
				const Eigen::Vector3d n0 = e.cross(e0); // / e.norm();
				const Eigen::Vector3d n1 = e1.cross(e); // / e.norm();
				orig_coeffs.block<1, 3>(k, 0) = vector_in_affine_coordinate(e, e0, n0, v);
				orig_coeffs.block<1, 3>(k, 3) = vector_in_affine_coordinate(e1, e, n1, v);
			}
		}
	}

	double SimilarityForm::value_unweighted(const Eigen::VectorXd &x) const
	{
		const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		double result = 0;
		for (int i = 0, k = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++, k++)
			{
				if (TT(i, j) < 0)
					continue;

				const Eigen::RowVectorXd &coeff = orig_coeffs.row(k);
				const Eigen::Vector3d e = (V.row(F_(i, le(j, 1))) - V.row(F_(i, le(j, 0))));
				const Eigen::Vector3d e0 = (V.row(F_(i, lv(j))) - V.row(F_(i, le(j, 0))));
				const Eigen::Vector3d e1 = (V.row(F_(TT(i, j), lv(TTi(i, j)))) - V.row(F_(i, le(j, 0))));
				const Eigen::Vector3d n0 = e.cross(e0); // / e.norm();
				const Eigen::Vector3d n1 = e1.cross(e); // / e.norm();
				const Eigen::Vector3d v0 = coeff(0) * e + coeff(1) * e0 + coeff(2) * n0;
				const Eigen::Vector3d v1 = coeff(3) * e1 + coeff(4) * e + coeff(5) * n1;

				const double err = (v0 - v1).squaredNorm();
				result += err * (orig_areas(i) + orig_areas(TT(i, j)));
			}
		}

		return result / 2.;
	}

	void SimilarityForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	{
		const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		gradv.setZero(V.size());
		for (int i = 0, k = 0; i < TT.rows(); i++)
		{
			for (int j = 0; j < TT.cols(); j++, k++)
			{
				if (TT(i, j) < 0)
					continue;

				Eigen::Matrix<double, 1, 18> tmp;
				tmp << V.row(F_(i, le(j, 0))), V.row(F_(i, le(j, 1))), V.row(F_(i, lv(j))), V.row(F_(TT(i, j), lv(TTi(i, j)))), orig_coeffs.row(k);
                Eigen::Vector<double, 12> g;
				similarity_gradient(tmp(0), tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11), tmp(12), tmp(13), tmp(14), tmp(15), tmp(16), tmp(17), g.data());

                g *= 0.5 * (orig_areas(i) + orig_areas(TT(i, j)));
                Eigen::Vector4i indices;
                indices << F_(i, le(j, 0)), F_(i, le(j, 1)), F_(i, lv(j)), F_(TT(i, j), lv(TTi(i, j)));
                for (int l = 0; l < 4; l++)
                    gradv.segment<3>(3 * indices(l)) += g.segment<3>(3 * l);
			}
		}
	}

	void SimilarityForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
	{
		POLYFEM_SCOPED_TIMER("similarity hessian");
		const Eigen::MatrixXd V = utils::unflatten(x, 3) + V_;

		Eigen::Matrix<int, 3, 2> le;
		le << 0, 1,
			1, 2,
			2, 0;
		Eigen::Vector3i lv;
		lv << 2, 0, 1;

		std::vector<Eigen::Triplet<double>> triplets;
		for (int i = 0, k = 0; i < TT.rows(); i++)
		{
			Eigen::Matrix<double, 1, 18> tmp;
			Eigen::Matrix<double, 12, 12> h;
			for (int j = 0; j < TT.cols(); j++, k++)
			{
				if (TT(i, j) < 0)
					continue;

				tmp << V.row(F_(i, le(j, 0))), V.row(F_(i, le(j, 1))), V.row(F_(i, lv(j))), V.row(F_(TT(i, j), lv(TTi(i, j)))), orig_coeffs.row(k);
				similarity_hessian(tmp(0), tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11), tmp(12), tmp(13), tmp(14), tmp(15), tmp(16), tmp(17), h.data());

				if (is_project_to_psd())
					h = ipc::project_to_psd(h);

                h *= 0.5 * (orig_areas(i) + orig_areas(TT(i, j)));
                Eigen::Vector4i indices;
                indices << F_(i, le(j, 0)), F_(i, le(j, 1)), F_(i, lv(j)), F_(TT(i, j), lv(TTi(i, j)));
                for (int lx = 0; lx < 4; lx++)
                    for (int ly = 0; ly < 4; ly++)
                        for (int dx = 0; dx < 3; dx++)
                            for (int dy = 0; dy < 3; dy++)
                                triplets.emplace_back(indices(lx) * 3 + dx, indices(ly) * 3 + dy, h(lx * 3 + dx, ly * 3 + dy));
			}
		}

		hessian.setZero();
		hessian.resize(x.size(), x.size());
		hessian.setFromTriplets(triplets.begin(), triplets.end());
	}
	// GlobalPositionalForm::GlobalPositionalForm(
	//     const Eigen::MatrixXd &V,
	//     const Eigen::MatrixXi &F,
	//     const Eigen::MatrixXd &source_skeleton_v,
	//     const Eigen::MatrixXd &target_skeleton_v,
	//     const Eigen::MatrixXi &skeleton_edges,
	//     const Eigen::MatrixXd &skin_weights):
	//     V_(V), source_skeleton_v_(source_skeleton_v),
	//     target_skeleton_v_(target_skeleton_v), skeleton_edges_(skeleton_edges),
	//     skin_weights_(skin_weights)
	// {
	//     Eigen::MatrixXd weights_per_bone = Eigen::MatrixXd::Zero(skeleton_edges.rows(), V.rows());
	//     for (int i = 0; i < skeleton_edges.rows(); i++)
	//         weights_per_bone.row(i) = skin_weights.row(skeleton_edges(i, 0)) + skin_weights.row(skeleton_edges(i, 1));

	//     bones = -Eigen::VectorXi::Ones(V.rows());
	//     relative_positions = Eigen::VectorXd::Zero(V.rows());
	//     for (int vid = 0; vid < V.rows(); vid++)
	//     {
	//         const Eigen::Vector3d point = V.row(vid);

	//         Eigen::MatrixXd::Index bid;
	//         if (weights_per_bone.col(vid).maxCoeff(&bid) < 0.8)
	//             continue;

	//         Eigen::Vector2d tmp = point_line_closest_distance(point, source_skeleton_v.row(skeleton_edges(bid, 0)), source_skeleton_v.row(skeleton_edges(bid, 1)));

	//         bones(vid) = bid;
	//         relative_positions(vid) = tmp(1);
	//     }
	// }

	// double GlobalPositionalForm::value_unweighted(const Eigen::VectorXd &x) const
	// {
	//     const double t = x(0);
	//     const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());
	//     const Eigen::MatrixXd skeleton_v = source_skeleton_v_ + t * (target_skeleton_v_ - source_skeleton_v_);

	//     double val = 0.;
	//     for (int vid = 0; vid < V.rows(); vid++)
	//     {
	//         const Eigen::Vector3d point = V.row(vid);

	//         if (bones(vid) < 0)
	//             continue;

	//         const int bid = bones(vid);
	//         const double param0 = relative_positions(vid);

	//         const Eigen::Vector2d tmp = point_line_closest_distance(point, skeleton_v.row(skeleton_edges_(bid, 0)), skeleton_v.row(skeleton_edges_(bid, 1)));
	//         val += pow(param0 - tmp(1), 2);
	//     }

	//     return val / 2.;
	// }

	// void GlobalPositionalForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	// {
	//     const double t = x(0);
	//     const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());
	//     const Eigen::MatrixXd skeleton_v = source_skeleton_v_ + t * (target_skeleton_v_ - source_skeleton_v_);

	//     gradv.setZero(x.size());
	//     for (int vid = 0; vid < V.rows(); vid++)
	//     {
	//         const Eigen::Vector3d point = V.row(vid);

	//         if (bones(vid) < 0)
	//             continue;

	//         const int bid = bones(vid);
	//         const double param0 = relative_positions(vid);

	//         const Eigen::Vector2d tmp = point_line_closest_distance(point, skeleton_v.row(skeleton_edges_(bid, 0)), skeleton_v.row(skeleton_edges_(bid, 1)));

	//         Eigen::Vector4d g;
	//         autogen::line_projection_uv_gradient(
	//             t, point(0), point(1), point(2), g.data(),
	//             source_skeleton_v_(skeleton_edges_(bid, 0), 0), source_skeleton_v_(skeleton_edges_(bid, 0), 1), source_skeleton_v_(skeleton_edges_(bid, 0), 2),
	//             target_skeleton_v_(skeleton_edges_(bid, 0), 0), target_skeleton_v_(skeleton_edges_(bid, 0), 1), target_skeleton_v_(skeleton_edges_(bid, 0), 2),
	//             source_skeleton_v_(skeleton_edges_(bid, 1), 0), source_skeleton_v_(skeleton_edges_(bid, 1), 1), source_skeleton_v_(skeleton_edges_(bid, 1), 2),
	//             target_skeleton_v_(skeleton_edges_(bid, 1), 0), target_skeleton_v_(skeleton_edges_(bid, 1), 1), target_skeleton_v_(skeleton_edges_(bid, 1), 2));

	//         g *= tmp(1) - param0;

	//         gradv(0) += g(0);
	//         gradv.template segment<3>(1 + vid * 3) += g.tail<3>();
	//     }
	// }

	// void GlobalPositionalForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
	// {
	//     hessian.resize(x.size(), x.size());
	//     hessian.setZero();

	//     const double t = x(0);
	//     const Eigen::MatrixXd V = V_ + utils::unflatten(x.tail(x.size() - 1), V_.cols());
	//     const Eigen::MatrixXd skeleton_v = source_skeleton_v_ + t * (target_skeleton_v_ - source_skeleton_v_);

	//     std::vector<Eigen::Triplet<double>> T;
	//     for (int vid = 0; vid < V.rows(); vid++)
	//     {
	//         const Eigen::Vector3d point = V.row(vid);

	//         if (bones(vid) < 0)
	//             continue;

	//         const int bid = bones(vid);
	//         const double param0 = relative_positions(vid);

	//         const Eigen::Vector2d tmp = point_line_closest_distance(point, skeleton_v.row(skeleton_edges_(bid, 0)), skeleton_v.row(skeleton_edges_(bid, 1)));

	//         Eigen::Vector4d g;
	//         autogen::line_projection_uv_gradient(
	//             t, point(0), point(1), point(2), g.data(),
	//             source_skeleton_v_(skeleton_edges_(bid, 0), 0), source_skeleton_v_(skeleton_edges_(bid, 0), 1), source_skeleton_v_(skeleton_edges_(bid, 0), 2),
	//             target_skeleton_v_(skeleton_edges_(bid, 0), 0), target_skeleton_v_(skeleton_edges_(bid, 0), 1), target_skeleton_v_(skeleton_edges_(bid, 0), 2),
	//             source_skeleton_v_(skeleton_edges_(bid, 1), 0), source_skeleton_v_(skeleton_edges_(bid, 1), 1), source_skeleton_v_(skeleton_edges_(bid, 1), 2),
	//             target_skeleton_v_(skeleton_edges_(bid, 1), 0), target_skeleton_v_(skeleton_edges_(bid, 1), 1), target_skeleton_v_(skeleton_edges_(bid, 1), 2));

	//         Eigen::Matrix4d h;
	//         autogen::line_projection_uv_hessian(
	//             t, point(0), point(1), point(2), h.data(),
	//             source_skeleton_v_(skeleton_edges_(bid, 0), 0), source_skeleton_v_(skeleton_edges_(bid, 0), 1), source_skeleton_v_(skeleton_edges_(bid, 0), 2),
	//             target_skeleton_v_(skeleton_edges_(bid, 0), 0), target_skeleton_v_(skeleton_edges_(bid, 0), 1), target_skeleton_v_(skeleton_edges_(bid, 0), 2),
	//             source_skeleton_v_(skeleton_edges_(bid, 1), 0), source_skeleton_v_(skeleton_edges_(bid, 1), 1), source_skeleton_v_(skeleton_edges_(bid, 1), 2),
	//             target_skeleton_v_(skeleton_edges_(bid, 1), 0), target_skeleton_v_(skeleton_edges_(bid, 1), 1), target_skeleton_v_(skeleton_edges_(bid, 1), 2));

	//         h = h.eval() * (tmp(1) - param0) + g * g.transpose();

	//         if (is_project_to_psd())
	//             h = ipc::project_to_psd(h);

	//     std::array<int, 4> index_map{{0, 1 + vid * 3, 2 + vid * 3, 3 + vid * 3}};
	//     for (int d0 = 0; d0 < 4; d0++)
	//         for (int d1 = 0; d1 < 4; d1++)
	//             T.emplace_back(index_map[d0], index_map[d1], h(d0, d1));
	// }

	//     hessian.setFromTriplets(T.begin(), T.end());
	// }
} // namespace polyfem::solver
