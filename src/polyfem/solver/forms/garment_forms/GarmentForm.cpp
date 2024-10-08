#include "GarmentForm.hpp"
#include <igl/triangle_triangle_adjacency.h>
#include <polyfem/utils/Logger.hpp>
#include <polyfem/utils/Timer.hpp>
#include <polyfem/autogen/auto_derivatives.hpp>
#include <finitediff.hpp>

using namespace polyfem::autogen;

namespace {
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
}

namespace polyfem::solver {
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
                
                const Eigen::Matrix<double, 9, 9> real_local_hess = local_grad * (local_grad.transpose() * barrier_hess) + local_hess * barrier_grad;
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        for (int di = 0; di < 3; di++)
                            for (int dj = 0; dj < 3; dj++)
                                triplets.emplace_back(3 * F_(e, i) + di, 3 * F_(e, j) + dj, real_local_hess(i * 3 + di, j * 3 + dj));
            }
        }

        hessian.setFromTriplets(triplets.begin(), triplets.end());
    }

    AngleForm::AngleForm(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) : V_(V), F_(F)
    {
        igl::triangle_triangle_adjacency(F_, TT, TTi);

        Eigen::MatrixXd normals = compute_normals(V, F_);

        areas = normals.rowwise().norm() / 2.;

        Eigen::Matrix<int, 3, 2> le;
        le << 0, 1, 
              1, 2, 
              2, 0;
        Eigen::Vector3i lv;
        lv << 2, 0, 1;

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

                result += (errA*errA + errB*errB) * (areas(i) + areas(TT(i, j)));
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

    double SimilarityForm::value_unweighted(const Eigen::VectorXd &x) const
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

    void SimilarityForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
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

    void SimilarityForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
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
                g({0,1,2,3,4,5,6,7,8}) += (1. / orig_dists(i, 2 * j + 0)) * local_grad[0];
                g({0,1,2,3,4,5,9,10,11}) -= (1. / orig_dists(i, 2 * j + 1)) * local_grad[1];
                Eigen::Matrix<double, 12, 12> h = (fac * g) * g.transpose();

                h({0,1,2,3,4,5,6,7,8}, {0,1,2,3,4,5,6,7,8})     += (fac * err / orig_dists(i, 2 * j + 0)) * local_hess[0];
                h({0,1,2,3,4,5,9,10,11}, {0,1,2,3,4,5,9,10,11}) -= (fac * err / orig_dists(i, 2 * j + 1)) * local_hess[1];

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
}