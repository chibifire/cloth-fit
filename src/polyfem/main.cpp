#include <filesystem>

#include <CLI/CLI.hpp>

#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/edges.h>

#include <ipc/ipc.hpp>

#include <polysolve/nonlinear/Solver.hpp>

#include <polyfem/State.hpp>
#include <polyfem/solver/forms/ContactForm.hpp>
#include <polyfem/solver/forms/GarmentForm.hpp>
#include <polyfem/solver/FullNLProblem.hpp>
#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/utils/JSONUtils.hpp>
#include <polyfem/utils/Logger.hpp>

using namespace polyfem;
using namespace solver;

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

int forward_simulation(const CLI::App &command_line,
					   const std::string &hdf5_file,
					   const std::string output_dir,
					   const unsigned max_threads,
					   const bool is_strict,
					   const bool fallback_solver,
					   const spdlog::level::level_enum &log_level,
					   json &in_args);

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

	std::string root = "/Users/zizhouhuang/Desktop/cloth-fit/python_clothing_deformer/scripts/initialization/";
	std::string avatar_mesh_path = root + "avatar.obj";
	std::string skinny_avatar_mesh_path = root + "projected_standard.obj";
	std::string garment_mesh_path = root + "garment.obj";

	const double scaling = 1e2;

	Eigen::MatrixXd avatar_v;
	Eigen::MatrixXi avatar_f;
	igl::read_triangle_mesh(avatar_mesh_path, avatar_v, avatar_f);
	avatar_v *= scaling;

	Eigen::MatrixXd skinny_avatar_v;
	Eigen::MatrixXi skinny_avatar_f;
	igl::read_triangle_mesh(skinny_avatar_mesh_path, skinny_avatar_v, skinny_avatar_f);
	skinny_avatar_v *= scaling;

	{
		Eigen::Vector3d center = (skinny_avatar_v.colwise().sum() - avatar_v.colwise().sum()) / avatar_v.rows();
		avatar_v.rowwise() += center.transpose();
	}

	Eigen::MatrixXd garment_v;
	Eigen::MatrixXi garment_f;
	igl::read_triangle_mesh(garment_mesh_path, garment_v, garment_f);
	garment_v *= scaling;

	// remove duplicate vertices in the garment
	{
		Eigen::VectorXi svi, svj;
		Eigen::MatrixXi sf;
		Eigen::MatrixXd sv;
		igl::remove_duplicate_vertices(garment_v, garment_f, 1e-4, sv, svi, svj, sf);
		std::swap(sv, garment_v);
		std::swap(sf, garment_f);
	}

	logger().info("avatar n_verts: {}, garment n_verts: {}, total n_verts: {}", avatar_v.rows(), garment_v.rows(), avatar_v.rows() + garment_v.rows());

	ipc::CollisionMesh collision_mesh;
	Eigen::MatrixXd collision_vertices(skinny_avatar_v.rows() + garment_v.rows(), garment_v.cols());
	Eigen::MatrixXi collision_triangles(skinny_avatar_f.rows() + garment_f.rows(), garment_f.cols());
	Eigen::MatrixXi collision_edges;
	{
		collision_vertices << skinny_avatar_v, garment_v;
		collision_triangles << skinny_avatar_f, garment_f.array() + skinny_avatar_v.rows();
		igl::edges(collision_triangles, collision_edges);
		collision_triangles = garment_f.array() + skinny_avatar_v.rows();
		collision_mesh = ipc::CollisionMesh(
			collision_vertices, collision_edges, collision_triangles);

		collision_triangles.resize(skinny_avatar_f.rows() + garment_f.rows(), garment_f.cols());
		collision_triangles << skinny_avatar_f, garment_f.array() + skinny_avatar_v.rows();

		const int n_avatar_verts = skinny_avatar_v.rows();
		collision_mesh.can_collide = [n_avatar_verts](size_t vi, size_t vj) {
			return vi >= n_avatar_verts || vj >= n_avatar_verts;
		};

		auto ids = ipc::my_has_intersections(collision_mesh, collision_vertices, ipc::BroadPhaseMethod::BVH);
		if (ids[0] >= 0)
		{
			io::OBJWriter::write(
				"intersection.obj", collision_vertices,
				collision_mesh.edges(), collision_mesh.faces());
			Eigen::MatrixXi edge(1, 2);
			edge << ids[0], ids[1];
			Eigen::MatrixXi face(1, 3);
			face << ids[2], ids[3], ids[4];
			io::OBJWriter::write(
				"intersecting_pair.obj", collision_vertices,
				edge, face);
			log_and_throw_error("Unable to solve, initial solution has intersections!");
		}
	}

	std::shared_ptr<ContactForm> contact_form;
	std::vector<std::shared_ptr<Form>> forms;
	{
		const double dhat = state.args["contact"]["dhat"];
		contact_form = std::make_shared<ContactForm>(collision_mesh, dhat, 1, false, false, false, false, state.args["solver"]["contact"]["CCD"]["broad_phase"], state.args["solver"]["contact"]["CCD"]["tolerance"], state.args["solver"]["contact"]["CCD"]["max_iterations"]);
		contact_form->set_weight(1);
		contact_form->set_barrier_stiffness(state.args["solver"]["contact"]["barrier_stiffness"]);
		contact_form->vis_collision_mesh_ = ipc::CollisionMesh(
				collision_vertices, collision_edges, collision_triangles);

		std::vector<int> indices(avatar_v.size());
		for (int i = 0; i < indices.size(); i++)
			indices[i] = i;
		auto penalty_form = std::make_shared<PenaltyForm>(utils::flatten(avatar_v - skinny_avatar_v), indices);
		penalty_form->set_weight(1e2);

		forms.push_back(penalty_form);
	}

	forms.push_back(contact_form);
	FullNLProblem nl_problem(forms);

	Eigen::VectorXd sol = utils::flatten(collision_vertices);
	sol.setZero();

	auto nl_solver = polysolve::nonlinear::Solver::create(in_args["solver"]["nonlinear"], in_args["solver"]["linear"], 1, logger());

	igl::write_triangle_mesh("initial.obj", contact_form->compute_displaced_surface(sol), collision_triangles);

	nl_problem.line_search_begin(sol, sol);
	if (!std::isfinite(nl_problem.value(sol))
		|| !nl_problem.is_step_valid(sol, sol)
		|| !nl_problem.is_step_collision_free(sol, sol))
		log_and_throw_error("Failed to apply boundary conditions!");

	nl_problem.init(sol);
	try
	{
		nl_solver->minimize(nl_problem, sol);
	}
	catch (const std::runtime_error &e)
	{
		throw e;
	}

	igl::write_triangle_mesh("result.obj", contact_form->compute_displaced_surface(sol), collision_triangles);

	return EXIT_SUCCESS;
}
