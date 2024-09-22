////////////////////////////////////////////////////////////////////////////////
#include <polyfem/solver/forms/garment_forms/GarmentForm.hpp>
#include <polyfem/solver/forms/garment_forms/CurveConstraintForm.hpp>

#include <finitediff.hpp>

#include <polyfem/State.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <igl/read_triangle_mesh.h>

#include <iostream>
#include <memory>
////////////////////////////////////////////////////////////////////////////////

using namespace polyfem;
using namespace polyfem::solver;
using namespace polyfem::time_integrator;
using namespace polyfem::assembler;

namespace
{
	std::shared_ptr<State> get_state(int dim, const std::string &material_type = "NeoHookean")
	{
		const std::string path = POLYFEM_DATA_DIR;

		json material = R"(
			{
				"type": "NeoHookean",
				"E": 20000,
				"nu": 0.3,
				"rho": 1000
			}
			)"_json;

		json in_args = R"(
		{
			"output": {
				"log": {
					"level": "warning"
				}
			}
		})"_json;
		in_args["materials"] = material;
    
        in_args["geometry"] = R"([{
            "transformation": {
                "scale": [0.1, 1, 1]
            },
            "n_refs": 1
        }])"_json;
        in_args["geometry"][0]["mesh"] = path + "/contact/meshes/3D/simple/bar/bar-6.msh";

		auto state = std::make_shared<State>();
		state->init(in_args, true);
		state->set_max_threads(1);

		return state;
	}
} // namespace

TEST_CASE("Garment forms derivatives", "[form][form_derivatives][garment]")
{
	const int dim = 3;
	const auto state_ptr = get_state(dim);

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
    igl::read_triangle_mesh("/Users/zizhouhuang/Desktop/cloth-fit/cpp_clothing_deformer/garment.obj", V, F);
	
	std::vector<std::unique_ptr<Form>> forms;
    forms.push_back(std::make_unique<CurveCurvatureForm>(V, F));
	forms.push_back(std::make_unique<AngleForm>(V, F));
	forms.push_back(std::make_unique<SimilarityForm>(V, F));

	static const int n_rand = 10;
	const double step = 1e-8;
	const double tol = 1e-6;


	for (auto &form : forms)
	{
		Eigen::VectorXd x = Eigen::VectorXd::Zero(V.size());

		form->init(x);
		form->init_lagging(x);

        Eigen::VectorXd grad;
        form->first_derivative(x, grad);

        REQUIRE(grad.norm() < 1e-12);
    }

    forms.push_back(std::make_unique<AreaForm>(V, F, 1));

	for (auto &form : forms)
	{
		Eigen::VectorXd x = Eigen::VectorXd::Zero(V.size());

        x.setRandom();
        x /= 100;

		form->init(x);
		form->init_lagging(x);

		for (int rand = 0; rand < n_rand; ++rand)
		{
			// Test gradient with finite differences
			{
				Eigen::VectorXd grad;
				form->first_derivative(x, grad);

				Eigen::VectorXd fgrad;
				fd::finite_gradient(
					x, [&form](const Eigen::VectorXd &x) -> double { return form->value(x); }, fgrad,
					fd::AccuracyOrder::SECOND, step);

				std::cout << std::setprecision(12) << "grad: " << grad.norm() << ", fd: " << fgrad.norm() << "\n";
				std::cout << std::setprecision(12) << "relative error: " << (grad - fgrad).norm() / fgrad.norm() << "\n";
				REQUIRE((grad - fgrad).norm() < tol * std::max(1.e-4, fgrad.norm()));
				// std::cout << grad.transpose() << "\n\n" << fgrad.transpose() << "\n\n";
			}

			// Test hessian with finite differences
			{
				StiffnessMatrix hess;
				form->second_derivative(x, hess);

				Eigen::MatrixXd fhess;
				fd::finite_jacobian(
					x,
					[&form](const Eigen::VectorXd &x) -> Eigen::VectorXd {
						Eigen::VectorXd grad;
						form->first_derivative(x, grad);
						return grad;
					},
					fhess,
					fd::AccuracyOrder::SECOND, step);

				std::cout << std::setprecision(12) << "hess: " << hess.norm() << ", fd: " << fhess.norm() << "\n";
				std::cout << std::setprecision(12) << "relative error: " << (hess - fhess).norm() / fhess.norm() << "\n";
				REQUIRE((hess - fhess).norm() < tol * std::max(1.e-4, fhess.norm()));
				// std::cout << (Eigen::MatrixXd)hess << "\n\n" << fhess << "\n\n";
			}
		}
	}
}
