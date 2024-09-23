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
#include <polyfem/solver/GarmentNLProblem.hpp>
#include <polyfem/solver/ALSolver.hpp>
#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/utils/JSONUtils.hpp>
#include <polyfem/utils/Logger.hpp>

#include <fstream>

using namespace polyfem;
using namespace solver;

template <int dim>
class Transformation {
public:
	/// @brief y = A * x + b
	/// @param A 
	/// @param b 
	Transformation(const Eigen::Matrix<double, dim, dim> &A, const Eigen::Vector<double, dim> &b) : A_(A), b_(b)
	{
	}

	void apply(Eigen::MatrixXd &V) const
	{
		V = (V * A_.transpose()).eval().rowwise() + b_.transpose();
	}

	void invert(Eigen::MatrixXd &V) const
	{
		V = (V.rowwise() - b_).eval() * A_.transpose().inv();
	}

private:
	const Eigen::Matrix<double, dim, dim> A_;
	const Eigen::Vector<double, dim> b_;
};

void read_edge_mesh(
	const std::string &path,
	Eigen::MatrixXd &V,
	Eigen::MatrixXi &E)
{
	std::ifstream infile(path);

	V.resize(0, 3);
	E.resize(0, 2);
	std::string line;
	while (std::getline(infile, line))
	{
		std::istringstream iss(line.substr(2));
		if (utils::StringUtils::startswith(line, "v "))
		{
			V.conservativeResize(V.rows() + 1, 3);

			double a, b, c;
			if (iss >> a >> b >> c)
				V.row(V.rows()-1) << a, b, c;
			else
				log_and_throw_error("read_edge_mesh failed to load vertex {}", V.rows() - 1);
		}
		else if (utils::StringUtils::startswith(line, "l "))
		{
			E.conservativeResize(E.rows() + 1, 2);

			int a, b;
			if (iss >> a >> b)
				E.row(E.rows()-1) << a - 1, b - 1;
			else
				log_and_throw_error("read_edge_mesh failed to load edge {}", E.rows() - 1);
		}
		else if (utils::StringUtils::startswith(line, "f "))
		{
			log_and_throw_error("read_edge_mesh does not support faces!");
		}
	}
}

Eigen::Vector2d point_edge_closest_distance(
	const Eigen::Vector3d &p,
	const Eigen::Vector3d &a,
	const Eigen::Vector3d &b)
{
	const Eigen::Vector3d e = b - a;
	const Eigen::Vector3d d = p - a;
	const double t = e.dot(d) / e.squaredNorm();
	const double dist = (d - t * e).squaredNorm();
	return Eigen::Vector2d(dist, t);
}

Eigen::MatrixXd extract_curve_center_targets(
	const Eigen::MatrixXd &garment_v,
	const std::vector<Eigen::VectorXi> &curves,
	const Eigen::MatrixXd &skeleton_v,
	const Eigen::MatrixXi &skeleton_bones,
	const Eigen::MatrixXd &target_skeleton_v)
{
	Eigen::MatrixXd targets(curves.size(), 3);
	for (int j = 0; j < curves.size(); j++)
	{
		// Compute centers of curves
		Eigen::Vector3d center = garment_v(curves[j], Eigen::all).colwise().sum() / curves[j].size();

		// Project centers to original skeleton bones
		int id = 0;
		double closest_dist = std::numeric_limits<double>::max(), closest_uv = 0;
		for (int i = 0; i < skeleton_bones.rows(); i++)
		{
			Eigen::Vector2d tmp = point_edge_closest_distance(center, skeleton_v.row(skeleton_bones(i, 0)), skeleton_v.row(skeleton_bones(i, 1)));
			if (tmp(0) < closest_dist)
			{
				closest_dist = tmp(0);
				closest_uv = tmp(1);
				id = i;
			}
		}

		// Map positions to new skeleton bones
		targets.row(j) = closest_uv * (target_skeleton_v.row(skeleton_bones(id, 1)) - target_skeleton_v.row(skeleton_bones(id, 0))) + target_skeleton_v.row(skeleton_bones(id, 0));
	}

	return targets;
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

	const std::string avatar_mesh_path = in_args["avatar_mesh_path"];
	const std::string skinny_avatar_mesh_path = in_args["skinny_avatar_mesh_path"];
	const std::string garment_mesh_path = in_args["garment_mesh_path"];

	const double scaling = 1e2;

	Eigen::MatrixXd avatar_v;
	Eigen::MatrixXi avatar_f;
	igl::read_triangle_mesh(avatar_mesh_path, avatar_v, avatar_f);

	Eigen::MatrixXd skinny_avatar_v;
	Eigen::MatrixXi skinny_avatar_f;
	igl::read_triangle_mesh(skinny_avatar_mesh_path, skinny_avatar_v, skinny_avatar_f);

	const Eigen::Vector3d center = (skinny_avatar_v.colwise().sum() - avatar_v.colwise().sum()) / avatar_v.rows();
	Transformation<3> trans(scaling * Eigen::Matrix3d::Identity(), scaling * center);

	skinny_avatar_v *= scaling;
	trans.apply(avatar_v);

	skinny_avatar_v += (avatar_v - skinny_avatar_v) * 1e-4;

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
	std::shared_ptr<PointPenaltyForm> pen_form;
	std::shared_ptr<PointLagrangianForm> lagr_form;
	std::vector<std::shared_ptr<Form>> forms;
	{
		const double dhat = state.args["contact"]["dhat"];
		contact_form = std::make_shared<ContactForm>(collision_mesh, dhat, 1, false, false, false, false, state.args["solver"]["contact"]["CCD"]["broad_phase"], state.args["solver"]["contact"]["CCD"]["tolerance"], state.args["solver"]["contact"]["CCD"]["max_iterations"]);
		contact_form->set_weight(1);
		contact_form->set_barrier_stiffness(state.args["solver"]["contact"]["barrier_stiffness"]);
		contact_form->save_ccd_debug_meshes = state.args["output"]["advanced"]["save_ccd_debug_meshes"];

		std::vector<int> indices(avatar_v.size());
		for (int i = 0; i < indices.size(); i++)
			indices[i] = i;
		pen_form = std::make_shared<PointPenaltyForm>(utils::flatten(avatar_v - skinny_avatar_v), indices);
		forms.push_back(pen_form);

		lagr_form = std::make_shared<PointLagrangianForm>(utils::flatten(avatar_v - skinny_avatar_v), indices);
		forms.push_back(lagr_form);

		auto angle_form = std::make_shared<AngleForm>(collision_vertices, collision_triangles.bottomRows(garment_f.rows()));
		angle_form->set_weight(state.args["angle_penalty_weight"]);
		forms.push_back(angle_form);

		auto similarity_form = std::make_shared<SimilarityForm>(collision_vertices, collision_triangles.bottomRows(garment_f.rows()));
		similarity_form->set_weight(state.args["similarity_penalty_weight"]);
		forms.push_back(similarity_form);

		auto curves = boundary_curves(collision_triangles.bottomRows(garment_f.rows()));

		auto curvature_form = std::make_shared<CurveCurvatureForm>(collision_vertices, curves);
		curvature_form->set_weight(state.args["curvature_penalty_weight"]);
		forms.push_back(curvature_form);

		{
			Eigen::MatrixXd skeleton_v, target_skeleton_v;
			Eigen::MatrixXi skeleton_bones, target_skeleton_bones;
			read_edge_mesh(state.args["original_skeleton_path"], skeleton_v, skeleton_bones);
			read_edge_mesh(state.args["target_skeleton_path"], target_skeleton_v, target_skeleton_bones);
			skeleton_v *= scaling;

			assert((skeleton_bones - target_skeleton_bones).squaredNorm() < 1);

			Eigen::MatrixXd centers = extract_curve_center_targets(collision_vertices, curves, skeleton_v, skeleton_bones, target_skeleton_v);
			trans.apply(centers);

			std::cout << centers << std::endl;
			
			auto center_target_form = std::make_shared<CurveCenterTargetForm>(collision_vertices, curves, centers);
			center_target_form->set_weight(state.args["curve_center_target_weight"]);
			forms.push_back(center_target_form);
		}
	}

	forms.push_back(contact_form);
	GarmentNLProblem nl_problem(1 + garment_v.size(), avatar_v - skinny_avatar_v, forms);

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
		[&](const Eigen::VectorXd &x) {
			state.solve_data.update_barrier_stiffness(sol);
		});

	// al_solver.post_subsolve = [&](const double al_weight) {
	// 	stats.solver_info.push_back(
	// 		{{"type", al_weight > 0 ? "al" : "rc"},
	// 			{"t", t}, // TODO: null if static?
	// 			{"info", nl_solver->info()}});
	// 	if (al_weight > 0)
	// 		stats.solver_info.back()["weight"] = al_weight;
	// 	save_subsolve(++subsolve_count, t, sol, Eigen::MatrixXd()); // no pressure
	// };

	Eigen::MatrixXd prev_sol = sol;
	al_solver.solve_al(nl_solver, nl_problem, sol);

	// igl::write_triangle_mesh("resultA.obj", contact_form->compute_displaced_surface(sol), collision_triangles);

	nl_solver = polysolve::nonlinear::Solver::create(state.args["solver"]["nonlinear"], state.args["solver"]["linear"], 1., logger());
	al_solver.solve_reduced(nl_solver, nl_problem, sol);

	// igl::write_triangle_mesh("resultB.obj", contact_form->compute_displaced_surface(sol), collision_triangles);

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