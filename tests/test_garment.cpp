////////////////////////////////////////////////////////////////////////////////
#include <polyfem/solver/forms/garment_forms/GarmentForm.hpp>
#include <polyfem/solver/forms/garment_forms/CurveConstraintForm.hpp>
#include <polyfem/solver/forms/garment_forms/CurveCenterTargetForm.hpp>
#include <polyfem/solver/forms/garment_forms/FitForm.hpp>

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
	
	auto curves = boundary_curves(F);
	Eigen::MatrixXd target(curves.size(), 3);
	target.setRandom();

	std::vector<std::unique_ptr<Form>> forms;
    forms.push_back(std::make_unique<CurveCurvatureForm>(V, curves));
	forms.push_back(std::make_unique<AngleForm>(V, F));
	forms.push_back(std::make_unique<SimilarityForm>(V, F));

	static const int n_rand = 10;
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

	forms.push_back(std::make_unique<CurveCenterTargetForm>(V, curves, target));
    forms.push_back(std::make_unique<AreaForm>(V, F, 1));
	
	{
		Eigen::MatrixXd avatar_v;
		Eigen::MatrixXi avatar_f;
		igl::read_triangle_mesh("/Users/zizhouhuang/Desktop/cloth-fit/cpp_clothing_deformer/avatar.obj", avatar_v, avatar_f);
		
		forms.push_back(std::make_unique<FitForm>(V, avatar_v, avatar_f, 0, 0.2));
	}

	for (auto &form : forms)
	{
		Eigen::VectorXd x = Eigen::VectorXd::Zero(V.size());

		for (int rand = 0; rand < n_rand; ++rand)
		{
			x.setRandom();
			x /= 100;

			double step;
			if (form->name() != "garment-fit")
				step = 1e-8;
			else
				step = 1e-6;

			form->init(x);
			form->init_lagging(x);
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
