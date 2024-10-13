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
	std::unordered_set<std::string> existing_names;
	for (const auto &form : prob.forms())
	{
		Eigen::VectorXd grad;
		form->first_derivative(complete_disp, grad);
		std::string name = "grad_" + form->name();
		while (existing_names.count(name) != 0)
			name += "_";
		existing_names.insert(name);
		writer.add_field(name, utils::unflatten(grad, 3));
		total_grad += grad;
	}
	writer.add_field("grad", utils::unflatten(total_grad, 3));

	Eigen::VectorXd body_ids = Eigen::VectorXd::Zero(V.rows());
	body_ids.head(n_avatar_vertices).array() = 1;
	writer.add_field("body_ids", body_ids);

	writer.write_mesh(path, utils::unflatten(complete_disp, V.cols()) + V, F);
}

Eigen::Vector3d bbox_size(const Eigen::Matrix<double, -1, 3> &V)
{
	return V.colwise().maxCoeff() - V.colwise().minCoeff();
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

	const std::string out_folder = state.args["/output/directory"_json_pointer];
	const std::string avatar_mesh_path = state.args["avatar_mesh_path"];
	const std::string skinny_avatar_mesh_path = state.args["skinny_avatar_mesh_path"];
	const std::string garment_mesh_path = state.args["garment_mesh_path"];
	const std::string source_skeleton_path = state.args["source_skeleton_path"];
	const std::string target_skeleton_path = state.args["target_skeleton_path"];

	// To scale the source mannequin and garment to around unit size
	const double source_scaling = 1e2;

	Eigen::MatrixXd avatar_v;
	Eigen::MatrixXi avatar_f;
	igl::read_triangle_mesh(avatar_mesh_path, avatar_v, avatar_f);

	Eigen::MatrixXd skeleton_v, target_skeleton_v;
	Eigen::MatrixXi skeleton_bones, target_skeleton_bones;
	read_edge_mesh(source_skeleton_path, skeleton_v, skeleton_bones);
	read_edge_mesh(target_skeleton_path, target_skeleton_v, target_skeleton_bones);
	skeleton_v *= source_scaling;

	assert((skeleton_bones - target_skeleton_bones).squaredNorm() < 1);

	Eigen::MatrixXd garment_v;
	Eigen::MatrixXi garment_f;
	{
		igl::read_triangle_mesh(garment_mesh_path, garment_v, garment_f);
		int n_refs = state.args["geometry"][0]["n_refs"];
		while (n_refs-- > 0)
			std::tie(garment_v, garment_f) = refine(garment_v, garment_f);
		garment_v *= source_scaling;

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

	Eigen::MatrixXd skinny_avatar_v;
	Eigen::MatrixXi skinny_avatar_f;
	igl::read_triangle_mesh(skinny_avatar_mesh_path, skinny_avatar_v, skinny_avatar_f);

	skinny_avatar_v *= source_scaling;

	const double target_scaling = bbox_size(skinny_avatar_v).maxCoeff() / bbox_size(target_skeleton_v).maxCoeff();
	const Eigen::Vector3d center = (skinny_avatar_v.colwise().sum() - target_scaling * avatar_v.colwise().sum()) / avatar_v.rows();
	Transformation<3> trans(target_scaling * Eigen::Matrix3d::Identity(), center);

	trans.apply(avatar_v);
	trans.apply(target_skeleton_v);

	skinny_avatar_v += (avatar_v - skinny_avatar_v) * 1e-4;

	igl::write_triangle_mesh(out_folder + "/target_avatar.obj", avatar_v, avatar_f);
	igl::write_triangle_mesh(out_folder + "/projected_avatar.obj", skinny_avatar_v, avatar_f);
	write_edge_mesh(out_folder + "/target_skeleton.obj", target_skeleton_v, target_skeleton_bones);
	write_edge_mesh(out_folder + "/source_skeleton.obj", skeleton_v, skeleton_bones);
	igl::write_triangle_mesh(out_folder + "/garment.obj", garment_v, garment_f);

	logger().info("avatar n_verts: {}, garment n_verts: {}, total n_verts: {}", avatar_v.rows(), garment_v.rows(), avatar_v.rows() + garment_v.rows());

	Eigen::MatrixXi collision_triangles(skinny_avatar_f.rows() + garment_f.rows(), garment_f.cols());
	collision_triangles << skinny_avatar_f, garment_f.array() + skinny_avatar_v.rows();
	Eigen::MatrixXi collision_edges;
	igl::edges(collision_triangles, collision_edges);

	Eigen::MatrixXd collision_vertices(skinny_avatar_v.rows() + garment_v.rows(), garment_v.cols());
	collision_vertices << skinny_avatar_v, garment_v;

	// std::cout << "target skeleton\n";
	// std::cout << target_skeleton_v << std::endl;

	// std::cout << "source skeleton\n";
	// std::cout << skeleton_v << std::endl;

	const auto curves = boundary_curves(collision_triangles.bottomRows(garment_f.rows()));
	const Eigen::MatrixXd source_curve_centers = extract_curve_center_targets(collision_vertices, curves, skeleton_v, skeleton_bones, skeleton_v);
	const Eigen::MatrixXd target_curve_centers = extract_curve_center_targets(collision_vertices, curves, skeleton_v, skeleton_bones, target_skeleton_v);

	// std::cout << "target curve centers\n";
	// std::cout << target_curve_centers << std::endl;

	// std::cout << "source curve centers\n";
	// std::cout << source_curve_centers << std::endl;

	Eigen::MatrixXd cur_garment_v = garment_v;
	int save_id = 1;
	const int total_steps = state.args["incremental_steps"];
	for (int substep = 0; substep < total_steps; ++substep)
	{
		const double prev_alpha = substep / (double)total_steps;
		const double next_alpha = (substep + 1) / (double)total_steps;

		// continuation
		const Eigen::MatrixXd prev_skeleton_v = (target_skeleton_v - skeleton_v) * prev_alpha + skeleton_v;
		const Eigen::MatrixXd next_skeleton_v = (target_skeleton_v - skeleton_v) * next_alpha + skeleton_v;
		const Eigen::MatrixXd prev_avatar_v = (avatar_v - skinny_avatar_v) * prev_alpha + skinny_avatar_v;
		const Eigen::MatrixXd next_avatar_v = (avatar_v - skinny_avatar_v) * next_alpha + skinny_avatar_v;

		const Eigen::MatrixXd next_curve_centers = (target_curve_centers - source_curve_centers) * next_alpha + source_curve_centers;

		ipc::CollisionMesh collision_mesh;
		{
			collision_vertices << prev_avatar_v, cur_garment_v;
			collision_mesh = ipc::CollisionMesh(
				collision_vertices, collision_edges, collision_triangles);

			const int n_avatar_verts = skinny_avatar_v.rows();
			collision_mesh.can_collide = [n_avatar_verts](size_t vi, size_t vj) {
				return vi >= n_avatar_verts || vj >= n_avatar_verts;
			};

			auto ids = ipc::my_has_intersections(collision_mesh, collision_vertices, ipc::BroadPhaseMethod::BVH);
			if (ids[0] >= 0)
			{
				io::OBJWriter::write(
					out_folder + "/intersection.obj", collision_vertices,
					collision_mesh.edges(), collision_mesh.faces());
				Eigen::MatrixXi edge(1, 2);
				edge << ids[0], ids[1];
				Eigen::MatrixXi face(1, 3);
				face << ids[2], ids[3], ids[4];
				io::OBJWriter::write(
					out_folder + "/intersecting_pair.obj", collision_vertices,
					edge, face);
				log_and_throw_error("Unable to solve, initial solution has intersections!");
			}
		}

		std::shared_ptr<ContactForm> contact_form;
		std::shared_ptr<PointPenaltyForm> pen_form;
		std::shared_ptr<PointLagrangianForm> lagr_form;
		std::shared_ptr<FitForm<4>> fit_form;
		std::vector<std::shared_ptr<Form>> forms;
		{
			const double dhat = state.args["contact"]["dhat"];
			contact_form = std::make_shared<ContactForm>(collision_mesh, dhat, 1, false, false, false, false, state.args["solver"]["contact"]["CCD"]["broad_phase"], state.args["solver"]["contact"]["CCD"]["tolerance"], state.args["solver"]["contact"]["CCD"]["max_iterations"]);
			contact_form->set_weight(1);
			contact_form->set_barrier_stiffness(state.args["solver"]["contact"]["barrier_stiffness"]);
			contact_form->save_ccd_debug_meshes = state.args["output"]["advanced"]["save_ccd_debug_meshes"];
			forms.push_back(contact_form);

			std::vector<int> indices(avatar_v.size());
			for (int i = 0; i < indices.size(); i++)
				indices[i] = i;
			pen_form = std::make_shared<PointPenaltyForm>(utils::flatten(next_avatar_v - prev_avatar_v), indices);
			forms.push_back(pen_form);

			lagr_form = std::make_shared<PointLagrangianForm>(utils::flatten(next_avatar_v - prev_avatar_v), indices);
			forms.push_back(lagr_form);

			auto angle_form = std::make_shared<AngleForm>(collision_vertices, collision_triangles.bottomRows(garment_f.rows()));
			angle_form->set_weight(state.args["angle_penalty_weight"]);
			forms.push_back(angle_form);

			auto similarity_form = std::make_shared<SimilarityForm>(collision_vertices, collision_triangles.bottomRows(garment_f.rows()));
			similarity_form->set_weight(state.args["similarity_penalty_weight"]);
			forms.push_back(similarity_form);

			auto curvature_form = std::make_shared<CurveCurvatureForm>(collision_vertices, curves);
			curvature_form->set_weight(state.args["curvature_penalty_weight"]);
			forms.push_back(curvature_form);

			// auto twist_form = std::make_shared<CurveTwistForm>(collision_vertices, curves);
			// twist_form->set_weight(state.args["twist_penalty_weight"]);
			// forms.push_back(twist_form);

			{
				auto center_target_form = std::make_shared<CurveCenterTargetForm>(collision_vertices, curves, next_curve_centers);
				center_target_form->set_weight(state.args["curve_center_target_weight"]);
				forms.push_back(center_target_form);
			}

			{
				fit_form = std::make_shared<FitForm<4>>(collision_vertices, collision_triangles.bottomRows(garment_f.rows()), avatar_v, avatar_f, 0.1);
				fit_form->disable();
				fit_form->set_weight(state.args["fit_weight"]);
				forms.push_back(fit_form);
			}

			for (int i = 0; i < curves.size(); i++)
			{
				auto form = std::make_shared<SymmetryForm>(collision_vertices, curves[i]);
				form->set_weight(state.args["symmetry_weight"]);
				if (form->enabled())
					forms.push_back(form);
			}
		}

		GarmentNLProblem nl_problem(1 + garment_v.size(), utils::flatten(next_avatar_v - prev_avatar_v), forms);

		Eigen::MatrixXd sol(nl_problem.full_size(), 1);
		sol.setZero();

		nl_problem.line_search_begin(sol, sol);
		if (!std::isfinite(nl_problem.value(sol))
			|| !nl_problem.is_step_valid(sol, sol)
			|| !nl_problem.is_step_collision_free(sol, sol))
			log_and_throw_error("Failed to apply boundary conditions!");


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
			const std::string path = out_folder + "/step_" + std::to_string(save_id++) + ".vtu";
			save_vtu(path, nl_problem, collision_vertices, collision_triangles, skinny_avatar_v.rows(), sol);
		};

		Eigen::MatrixXd prev_sol = sol;
		al_solver.solve_al(nl_solver, nl_problem, sol);

		fit_form->enable();

		nl_solver = polysolve::nonlinear::Solver::create(state.args["solver"]["nonlinear"], state.args["solver"]["linear"], 1., logger());
		al_solver.solve_reduced(nl_solver, nl_problem, sol);

		cur_garment_v += utils::unflatten(sol.bottomRows(cur_garment_v.size()), 3);
	}

	return EXIT_SUCCESS;
}

// #include <Eigen/core>
// #include <array>
// #include <iostream>

// template <int N>
// constexpr Eigen::Matrix<double, ((N+1)*(N+2))/2, 4> upsample_standard()
// {
// 	constexpr int num = ((N+1)*(N+2))/2;

// 	Eigen::Matrix<double, num, 4> out;
// 	for (int i = 0, k = 0; i <= N; i++)
// 		for (int j = 0; i + j <= N; j++, k++)
// 		{
// 			std::array<int, 3> arr = {i, j, N - i - j};
// 			std::sort(arr.begin(), arr.end());

// 			double w = 6;
// 			if (arr[1] == 0)        // vertex node
// 				w = 1;
// 			else if (arr[0] == 0)   // edge node
// 				w = 3;
// 			else                    // face node
// 				w = 6;
			
// 			out.row(k) << i, j, N - i - j, w;
// 		}
	
// 	out.template leftCols<3>() /= N;
// 	out.col(3) /= N * N;
// 	return out;
// }

// int main()
// {
// 	std::cout << upsample_standard<3>() << std::endl;
// }

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