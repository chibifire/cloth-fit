#include <filesystem>

#include <CLI/CLI.hpp>

#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/edges.h>

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
#include <polyfem/utils/JSONUtils.hpp>
#include <polyfem/utils/Logger.hpp>
#include <polyfem/mesh/MeshUtils.hpp>
#include <polyfem/garment/optimize.hpp>

#include <fstream>

using namespace polyfem;
using namespace solver;
using namespace mesh;

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

	GarmentSolver gstate;

	const std::string out_folder = state.args["/output/directory"_json_pointer];
	const std::string avatar_mesh_path = state.args["avatar_mesh_path"];
	const std::string garment_mesh_path = state.args["garment_mesh_path"];
	const std::string source_skeleton_path = state.args["source_skeleton_path"];
	const std::string target_skeleton_path = state.args["target_skeleton_path"];
	const std::string skin_weights_path = state.args["skin_weights_path"];

	gstate.out_folder = out_folder;

	gstate.read_meshes(avatar_mesh_path, source_skeleton_path, target_skeleton_path, skin_weights_path);
	gstate.load_garment_mesh(garment_mesh_path, state.args["geometry"][0]["n_refs"]);
	gstate.normalize_meshes();
	gstate.project_avatar_to_skeleton();

	igl::write_triangle_mesh(out_folder + "/target_avatar.obj", gstate.avatar_v, gstate.avatar_f);
	igl::write_triangle_mesh(out_folder + "/projected_avatar.obj", gstate.skinny_avatar_v, gstate.avatar_f);
	write_edge_mesh(out_folder + "/target_skeleton.obj", gstate.target_skeleton_v, gstate.target_skeleton_b);
	write_edge_mesh(out_folder + "/source_skeleton.obj", gstate.skeleton_v, gstate.skeleton_b);
	igl::write_triangle_mesh(out_folder + "/garment.obj", gstate.garment_v, gstate.garment_f);

	logger().info("avatar n_verts: {}, garment n_verts: {}, total n_verts: {}", gstate.avatar_v.rows(), gstate.garment_v.rows(), gstate.avatar_v.rows() + gstate.garment_v.rows());

	Eigen::MatrixXi collision_triangles(gstate.avatar_f.rows() + gstate.garment_f.rows(), gstate.garment_f.cols());
	collision_triangles << gstate.avatar_f, gstate.garment_f.array() + gstate.skinny_avatar_v.rows();
	Eigen::MatrixXi collision_edges;
	igl::edges(collision_triangles, collision_edges);

	Eigen::MatrixXd collision_vertices(gstate.skinny_avatar_v.rows() + gstate.garment_v.rows(), gstate.garment_v.cols());
	collision_vertices << gstate.skinny_avatar_v, gstate.garment_v;

	ipc::CollisionMesh collision_mesh;
	{
		collision_mesh = ipc::CollisionMesh(
			collision_vertices, collision_edges, collision_triangles);

		const int n_avatar_verts = gstate.skinny_avatar_v.rows();
		collision_mesh.can_collide = [n_avatar_verts](size_t vi, size_t vj) {
			return vi >= n_avatar_verts || vj >= n_avatar_verts;
		};

		gstate.check_intersections(collision_mesh, collision_vertices);
	}

	const auto curves = boundary_curves(collision_triangles.bottomRows(gstate.garment_f.rows()));
	const Eigen::MatrixXd source_curve_centers = extract_curve_center_targets(collision_vertices, curves, gstate.skeleton_v, gstate.skeleton_b, gstate.skeleton_v);
	const Eigen::MatrixXd target_curve_centers = extract_curve_center_targets(collision_vertices, curves, gstate.skeleton_v, gstate.skeleton_b, gstate.target_skeleton_v);

	Eigen::MatrixXd cur_garment_v = gstate.garment_v;
	int save_id = 0;
	const int total_steps = state.args["incremental_steps"];
	const int stride = state.args["output"]["paraview"]["skip_frame"];

	std::vector<std::shared_ptr<Form>> persistent_forms;
	std::shared_ptr<CurveSizeForm> curve_size_form;
	std::shared_ptr<ContactForm> contact_form;
	{
		auto angle_form = std::make_shared<AngleForm>(collision_vertices, collision_triangles.bottomRows(gstate.garment_f.rows()));
		angle_form->set_weight(state.args["angle_penalty_weight"]);
		persistent_forms.push_back(angle_form);

		auto def_form = std::make_shared<DefGradForm>(collision_vertices, collision_triangles.bottomRows(gstate.garment_f.rows()));
		def_form->set_weight(state.args["deformation_penalty_weight"]);
		persistent_forms.push_back(def_form);

		auto similarity_form = std::make_shared<SimilarityForm>(collision_vertices, collision_triangles.bottomRows(gstate.garment_f.rows()));
		similarity_form->set_weight(state.args["similarity_penalty_weight"]);
		persistent_forms.push_back(similarity_form);

		auto curvature_form = std::make_shared<CurveCurvatureForm>(collision_vertices, curves);
		curvature_form->set_weight(state.args["curvature_penalty_weight"]);
		persistent_forms.push_back(curvature_form);

		auto twist_form = std::make_shared<CurveTwistForm>(collision_vertices, curves);
		twist_form->set_weight(state.args["twist_penalty_weight"]);
		persistent_forms.push_back(twist_form);

		auto sym_form = std::make_shared<SymmetryForm>(collision_vertices, curves);
		sym_form->set_weight(state.args["symmetry_weight"]);
		if (sym_form->enabled())
			persistent_forms.push_back(sym_form);

		curve_size_form = std::make_shared<CurveSizeForm>(collision_vertices, curves);
		curve_size_form->disable();
		curve_size_form->set_weight(state.args["curve_size_weight"]);
		persistent_forms.push_back(curve_size_form);

		const double dhat = state.args["contact"]["dhat"];
		contact_form = std::make_shared<ContactForm>(collision_mesh, dhat, 1, false, false, false, false, state.args["solver"]["contact"]["CCD"]["broad_phase"], state.args["solver"]["contact"]["CCD"]["tolerance"], state.args["solver"]["contact"]["CCD"]["max_iterations"]);
		contact_form->set_weight(1);
		contact_form->set_barrier_stiffness(state.args["solver"]["contact"]["barrier_stiffness"]);
		contact_form->save_ccd_debug_meshes = state.args["output"]["advanced"]["save_ccd_debug_meshes"];
		persistent_forms.push_back(contact_form);
	}

	Eigen::MatrixXd sol = Eigen::MatrixXd::Zero(1 + gstate.garment_v.size(), 1);
	for (int substep = 0; substep < total_steps; ++substep)
	{
		const double prev_alpha = substep / (double)total_steps;
		const double next_alpha = (substep + 1) / (double)total_steps;

		logger().info("Start substep {} out of {}", substep + 1, total_steps);

		// continuation
		// const Eigen::MatrixXd prev_skeleton_v = (gstate.target_skeleton_v - gstate.skeleton_v) * prev_alpha + gstate.skeleton_v;
		// const Eigen::MatrixXd next_skeleton_v = (gstate.target_skeleton_v - gstate.skeleton_v) * next_alpha + gstate.skeleton_v;
		// const Eigen::MatrixXd prev_avatar_v = (gstate.avatar_v - gstate.skinny_avatar_v) * prev_alpha + gstate.skinny_avatar_v;
		const Eigen::MatrixXd next_avatar_v = (gstate.avatar_v - gstate.skinny_avatar_v) * next_alpha + gstate.skinny_avatar_v;

		const Eigen::MatrixXd next_curve_centers = (target_curve_centers - source_curve_centers) * next_alpha + source_curve_centers;

		std::vector<std::shared_ptr<Form>> forms = persistent_forms;
		std::shared_ptr<PointPenaltyForm> pen_form;
		std::shared_ptr<PointLagrangianForm> lagr_form;
		std::shared_ptr<FitForm<4>> fit_form;
		{
			std::vector<int> indices(gstate.avatar_v.size());
			for (int i = 0; i < indices.size(); i++)
				indices[i] = i;
			pen_form = std::make_shared<PointPenaltyForm>(utils::flatten(next_avatar_v - gstate.skinny_avatar_v), indices);
			forms.push_back(pen_form);

			lagr_form = std::make_shared<PointLagrangianForm>(utils::flatten(next_avatar_v - gstate.skinny_avatar_v), indices);
			forms.push_back(lagr_form);

			auto center_target_form = std::make_shared<CurveCenterTargetForm>(collision_vertices, curves, next_curve_centers);
			center_target_form->set_weight(state.args["curve_center_target_weight"]);
			forms.push_back(center_target_form);

			fit_form = std::make_shared<FitForm<4>>(collision_vertices, collision_triangles.bottomRows(gstate.garment_f.rows()), gstate.avatar_v, gstate.avatar_f, 0.1);
			fit_form->disable();
			fit_form->set_weight(state.args["fit_weight"]);
			forms.push_back(fit_form);

			curve_size_form->disable();
		}

		GarmentNLProblem nl_problem(1 + gstate.garment_v.size(), utils::flatten(gstate.avatar_v - gstate.skinny_avatar_v), forms, {});
		nl_problem.set_target_value(next_alpha);

		nl_problem.line_search_begin(sol, sol);
		if (!std::isfinite(nl_problem.value(sol))
			|| !nl_problem.is_step_valid(sol, sol)
			|| !nl_problem.is_step_collision_free(sol, sol))
			log_and_throw_error("Failed to apply boundary conditions!");

		std::shared_ptr<polysolve::nonlinear::Solver> nl_solver = polysolve::nonlinear::Solver::create(state.args["solver"]["augmented_lagrangian"]["nonlinear"], state.args["solver"]["linear"], 1., logger());

		double initial_weight = state.args["solver"]["augmented_lagrangian"]["initial_weight"];
		const double scaling = state.args["solver"]["augmented_lagrangian"]["scaling"];
		const double max_weight = state.args["solver"]["augmented_lagrangian"]["max_weight"].get<double>();
		while (true) {
			pen_form->set_weight(initial_weight);
			Eigen::VectorXd g;
			nl_problem.gradient(sol, g);
			if (g(0) > 0)
				break;
			initial_weight *= scaling;
			if (initial_weight >= max_weight)
			{
				initial_weight = max_weight;
				logger().debug("Increase initial AL weight to {}", initial_weight);
				break;
			}
		}

		ALSolver<GarmentNLProblem, PointLagrangianForm, PointPenaltyForm> al_solver(
			lagr_form, pen_form,
			initial_weight, scaling, max_weight,
			state.args["solver"]["augmented_lagrangian"]["eta"],
			state.args["solver"]["augmented_lagrangian"]["error_threshold"],
			[&](const Eigen::VectorXd &x) {
				state.solve_data.update_barrier_stiffness(sol);
			});

		nl_problem.post_step_call_back = [&](const Eigen::VectorXd &sol) {
			const std::string path = out_folder + "/step_" + std::to_string(save_id / stride) + ".vtu";
			if (save_id % stride == 0)
				save_vtu(path, nl_problem, collision_vertices, collision_triangles, gstate.skinny_avatar_v.rows(), sol);
			++save_id;
		};

		al_solver.solve_al(nl_solver, nl_problem, sol);

		fit_form->enable();
		if (substep == total_steps - 1)
			curve_size_form->enable();

		nl_solver = polysolve::nonlinear::Solver::create(state.args["solver"]["nonlinear"], state.args["solver"]["linear"], 1., logger());
		al_solver.solve_reduced(nl_solver, nl_problem, sol);

		cur_garment_v = gstate.garment_v + utils::unflatten(sol.bottomRows(cur_garment_v.size()), 3);
	}

	return EXIT_SUCCESS;
}
