#include <filesystem>

#include <CLI/CLI.hpp>

#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/edges.h>

#include <ipc/ipc.hpp>

#include <polysolve/nonlinear/Solver.hpp>

#include <polyfem/State.hpp>
#include <polyfem/utils/StringUtils.hpp>
#include <polyfem/solver/forms/ContactForm.hpp>
#include <polyfem/solver/forms/garment_forms/GarmentForm.hpp>
#include <polyfem/solver/forms/garment_forms/GarmentALForm.hpp>
#include <polyfem/solver/forms/garment_forms/CurveConstraintForm.hpp>
#include <polyfem/solver/forms/garment_forms/CurveCenterTargetForm.hpp>
#include <polyfem/solver/forms/garment_forms/FitForm.hpp>
#include <polyfem/solver/GarmentNLProblem.hpp>
#include <polyfem/solver/ALSolver.hpp>
#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/utils/JSONUtils.hpp>
#include <polyfem/utils/Logger.hpp>
#include <polyfem/mesh/MeshUtils.hpp>

#include <paraviewo/ParaviewWriter.hpp>
#include <paraviewo/VTUWriter.hpp>

#include <fstream>

using namespace polyfem;
using namespace solver;
using namespace mesh;

void save_vtu(
	const std::string &path,
	GarmentNLProblem &prob,
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const int n_avatar_vertices,
	const Eigen::VectorXd &sol)
{
	std::shared_ptr<paraviewo::ParaviewWriter> tmpw = std::make_shared<paraviewo::VTUWriter>();
	paraviewo::ParaviewWriter &writer = *tmpw;

	const Eigen::VectorXd complete_disp = prob.full_to_complete(prob.reduced_to_full(sol));

	Eigen::VectorXd total_grad = Eigen::VectorXd::Zero(complete_disp.size());
	for (const auto &form : prob.forms())
	{
		Eigen::VectorXd grad;
		form->first_derivative(complete_disp, grad);
		writer.add_field("grad_" + form->name(), utils::unflatten(grad, 3));
		total_grad += grad;
	}
	writer.add_field("grad", utils::unflatten(total_grad, 3));

	Eigen::VectorXd body_ids = Eigen::VectorXd::Zero(V.rows());
	body_ids.head(n_avatar_vertices).array() = 1;
	writer.add_field("body_ids", body_ids);

	writer.write_mesh(path, utils::unflatten(complete_disp, V.cols()) + V, F);
}

bool has_arg(const CLI::App &command_line, const std::string &value)
{
	const auto *opt = command_line.get_option_no_throw(value.size() == 1 ? ("-" + value) : ("--" + value));
	if (!opt)
		return false;

	return opt->count() > 0;
}

bool load_json(const std::string &json_file, json &out)
{
	std::ifstream file(json_file);

	if (!file.is_open())
		return false;

	file >> out;

	if (!out.contains("root_path"))
		out["root_path"] = json_file;

	return true;
}

int main(int argc, char **argv)
{
	using namespace polyfem;

	CLI::App command_line{"polyfem"};

	command_line.ignore_case();
	command_line.ignore_underscore();

	// Eigen::setNbThreads(1);
	unsigned max_threads = std::numeric_limits<unsigned>::max();
	command_line.add_option("--max_threads", max_threads, "Maximum number of threads");

	auto input = command_line.add_option_group("input");

	std::string json_file = "";
	input->add_option("-j,--json", json_file, "Simulation JSON file")->check(CLI::ExistingFile);

	input->require_option(1);

	std::string output_dir = "";
	command_line.add_option("-o,--output_dir", output_dir, "Directory for output files")->check(CLI::ExistingDirectory | CLI::NonexistentPath);

	const std::vector<std::pair<std::string, spdlog::level::level_enum>>
		SPDLOG_LEVEL_NAMES_TO_LEVELS = {
			{"trace", spdlog::level::trace},
			{"debug", spdlog::level::debug},
			{"info", spdlog::level::info},
			{"warning", spdlog::level::warn},
			{"error", spdlog::level::err},
			{"critical", spdlog::level::critical},
			{"off", spdlog::level::off}};
	spdlog::level::level_enum log_level = spdlog::level::debug;
	command_line.add_option("--log_level", log_level, "Log level")
		->transform(CLI::CheckedTransformer(SPDLOG_LEVEL_NAMES_TO_LEVELS, CLI::ignore_case));

	CLI11_PARSE(command_line, argc, argv);

	json in_args = json({});
	{
		if (!json_file.empty())
		{
			const bool ok = load_json(json_file, in_args);

			if (!ok)
				log_and_throw_error(fmt::format("unable to open {} file", json_file));
		}

		if (in_args.empty())
		{
			logger().error("No input file specified!");
			return command_line.exit(CLI::RequiredError("--json"));
		}

		json tmp = json::object();
		if (has_arg(command_line, "log_level"))
			tmp["/output/log/level"_json_pointer] = int(log_level);
		if (has_arg(command_line, "max_threads"))
			tmp["/solver/max_threads"_json_pointer] = max_threads;
		if (has_arg(command_line, "output_dir"))
			tmp["/output/directory"_json_pointer] = std::filesystem::absolute(output_dir);
		
		assert(tmp.is_object());
		in_args.merge_patch(tmp);
	}

	State state;
	state.init(in_args, false);

	const std::string garment_mesh_path = state.args["garment_mesh_path"];

	Eigen::MatrixXd garment_v;
	Eigen::MatrixXi garment_f;
	{
		igl::read_triangle_mesh(garment_mesh_path, garment_v, garment_f);
		int n_refs = state.args["geometry"][0]["n_refs"];
		while (n_refs-- > 0)
			std::tie(garment_v, garment_f) = refine(garment_v, garment_f);

		// remove duplicate vertices in the garment
		{
			Eigen::VectorXi svi, svj;
			Eigen::MatrixXi sf;
			Eigen::MatrixXd sv;
			igl::remove_duplicate_vertices(garment_v, garment_f, 1e-4, sv, svi, svj, sf);
			std::swap(sv, garment_v);
			std::swap(sf, garment_f);
		}
	}

	auto curves = boundary_curves(garment_f);

	Eigen::Vector3d min_, max_;
	min_ = garment_v.colwise().minCoeff();
	max_ = garment_v.colwise().maxCoeff();

	// reorder indices so that vertices on two ends are at the front
	std::vector<bool> mask(garment_v.rows(), false);
	{
		std::vector<int> indices;
		for (auto &c : curves)
		{
			for (auto i : c)
			{
				if (garment_v(i, 0) < min_(0) + 1e-5 || garment_v(i, 0) > max_(0) - 1e-5)
				{
					if (!mask[i])
						indices.push_back(i);
					mask[i] = true;
				}
			}
		}
		
		for (int i = 0; i < garment_v.rows(); i++)
		{
			if (!mask[i])
				indices.push_back(i);
		}

		garment_v = garment_v(indices, Eigen::all).eval();

		std::vector<int> inv_indices(indices.size());
		for (int i = 0; i < indices.size(); i++)
			inv_indices[indices[i]] = i;

		for (int f = 0; f < garment_f.rows(); f++)
		{
			for (int d = 0; d < 3; d++)
			{
				garment_f(f, d) = inv_indices[garment_f(f, d)];
				assert(garment_f(f, d) < garment_v.rows());
			}
		}

		std::vector<bool> new_mask = mask;
		for (int i = 0; i < mask.size(); i++)
		{
			new_mask[i] = mask[indices[i]];
		}
		std::swap(new_mask, mask);
	}

	ipc::CollisionMesh collision_mesh;
	{
		Eigen::MatrixXi collision_edges;
		igl::edges(garment_f, collision_edges);
		collision_mesh = ipc::CollisionMesh(
			garment_v, collision_edges, garment_f);

		auto ids = ipc::my_has_intersections(collision_mesh, garment_v, ipc::BroadPhaseMethod::BVH);
		if (ids[0] >= 0)
		{
			io::OBJWriter::write(
				"intersection.obj", garment_v,
				collision_mesh.edges(), collision_mesh.faces());
			Eigen::MatrixXi edge(1, 2);
			edge << ids[0], ids[1];
			Eigen::MatrixXi face(1, 3);
			face << ids[2], ids[3], ids[4];
			io::OBJWriter::write(
				"intersecting_pair.obj", garment_v,
				edge, face);
			log_and_throw_error("Unable to solve, initial solution has intersections!");
		}
	}

	Eigen::MatrixXd target_v = garment_v * 0.5;
	int save_id = 0;
	int stride = state.args["output"]["paraview"]["skip_frame"];
	const int total_steps = state.args["incremental_steps"];
	Eigen::MatrixXd sol;
	for (int substep = 0; substep < total_steps; ++substep)
	{
		const double next_alpha = (substep + 1) / (double)total_steps;

		std::vector<std::shared_ptr<Form>> forms;

		std::vector<int> indices;
		for (int i = 0; i < garment_v.rows(); i++)
			if (mask[i])
				for (int d = 0; d < 3; d++)
					indices.push_back(i * 3 + d);
		
		Eigen::VectorXd target_diff = next_alpha * utils::flatten(target_v - garment_v)(indices);

		auto pen_form = std::make_shared<PointPenaltyForm>(target_diff, indices);
		forms.push_back(pen_form);

		auto lagr_form = std::make_shared<PointLagrangianForm>(target_diff, indices);
		forms.push_back(lagr_form);

		auto angle_form = std::make_shared<AngleForm>(garment_v, garment_f);
		angle_form->set_weight(state.args["angle_penalty_weight"]);
		forms.push_back(angle_form);

		auto similarity_form = std::make_shared<SimilarityForm>(garment_v, garment_f);
		similarity_form->set_weight(state.args["similarity_penalty_weight"]);
		forms.push_back(similarity_form);

		auto curvature_form = std::make_shared<CurveCurvatureForm>(garment_v, curves);
		curvature_form->set_weight(state.args["curvature_penalty_weight"]);
		forms.push_back(curvature_form);

		// auto twist_form = std::make_shared<CurveTwistForm>(garment_v, curves);
		// twist_form->set_weight(state.args["twist_penalty_weight"]);
		// forms.push_back(twist_form);

		assert(true);

		if (state.args["contact"]["enabled"])
		{
			const double dhat = state.args["contact"]["dhat"];
			auto contact_form = std::make_shared<ContactForm>(collision_mesh, dhat, 1, false, false, false, false, state.args["solver"]["contact"]["CCD"]["broad_phase"], state.args["solver"]["contact"]["CCD"]["tolerance"], state.args["solver"]["contact"]["CCD"]["max_iterations"]);
			contact_form->set_weight(1);
			contact_form->set_barrier_stiffness(state.args["solver"]["contact"]["barrier_stiffness"]);
			contact_form->save_ccd_debug_meshes = state.args["output"]["advanced"]["save_ccd_debug_meshes"];
			forms.push_back(contact_form);
		}

		GarmentNLProblem nl_problem(1 + garment_v.size() - target_diff.size(), target_diff, forms);

		logger().info("Total nodes {}, number of constrained nodes {}", garment_v.size() / 3, target_diff.size() / 3);

		if (sol.size() == 0)
			sol.setZero(nl_problem.full_size(), 1);

		nl_problem.line_search_begin(sol, sol);
		if (!std::isfinite(nl_problem.value(sol))
			|| !nl_problem.is_step_valid(sol, sol)
			|| !nl_problem.is_step_collision_free(sol, sol))
			log_and_throw_error("Failed to apply boundary conditions!");

		nl_problem.line_search_begin(sol, sol);

		std::shared_ptr<polysolve::nonlinear::Solver> nl_solver = polysolve::nonlinear::Solver::create(state.args["solver"]["augmented_lagrangian"]["nonlinear"], state.args["solver"]["linear"], 1., logger());

		ALSolver<GarmentNLProblem, PointLagrangianForm, PointPenaltyForm> al_solver(
			lagr_form, pen_form,
			state.args["solver"]["augmented_lagrangian"]["initial_weight"],
			state.args["solver"]["augmented_lagrangian"]["scaling"],
			state.args["solver"]["augmented_lagrangian"]["max_weight"],
			state.args["solver"]["augmented_lagrangian"]["eta"],
			state.args["solver"]["augmented_lagrangian"]["error_threshold"],
			[&](const Eigen::VectorXd &x) {
				state.solve_data.update_barrier_stiffness(sol);
			});
		
		nl_problem.post_step_call_back = [&](const Eigen::VectorXd &sol) {
			if (save_id % stride == 0)
			{
				const std::string path = "step_" + std::to_string(save_id / stride) + ".vtu";
				save_vtu(path, nl_problem, garment_v, garment_f, target_diff.size() / 3, sol);
			}
			save_id++;
		};

		save_vtu("step_" + std::to_string(save_id++) + ".vtu", nl_problem, garment_v, garment_f, 0, sol);

		Eigen::MatrixXd prev_sol = sol;
		al_solver.solve_al(nl_solver, nl_problem, sol);

		nl_solver = polysolve::nonlinear::Solver::create(state.args["solver"]["nonlinear"], state.args["solver"]["linear"], 1., logger());
		al_solver.solve_reduced(nl_solver, nl_problem, sol);
	}

	return EXIT_SUCCESS;
}



// #include <iostream>
// #include <Eigen/SparseCore>

// void assign(Eigen::SparseMatrix<double, Eigen::ColMajor> &mat)
// {
// 	mat.resize(0, 0);

// 	Eigen::SparseMatrix<double, Eigen::ColMajor> A;
	
// 	A.resize(2206, 2206);
// 	A.setIdentity();

// 	mat = A;

// 	std::cout << "mat:\n" << mat.rows() << std::endl;
// 	std::cout << "A:\n" << A.rows() << std::endl;
// }

// int main() {
// 	Eigen::SparseMatrix<double, Eigen::ColMajor> B;
// 	assign(B);
// }