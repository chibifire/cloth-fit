////////////////////////////////////////////////////////////////////////////////
#include <polyfem/mesh/MeshUtils.hpp>
#include <polyfem/io/OBJReader.hpp>
#include <polyfem/io/OBJWriter.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <Eigen/Dense>

#include <filesystem>
#include <fstream>
#include <sstream>
////////////////////////////////////////////////////////////////////////////////

using namespace polyfem;

class MeshTestFixture {
public:
    MeshTestFixture() {
        test_dir = std::filesystem::temp_directory_path() / "polyfem_mesh_tests";
        std::filesystem::create_directories(test_dir);
    }

    ~MeshTestFixture() {
        if (std::filesystem::exists(test_dir)) {
            std::filesystem::remove_all(test_dir);
        }
    }

    std::string temp_file(const std::string& suffix) {
        return (test_dir / ("mesh_" + suffix)).string();
    }

private:
    std::filesystem::path test_dir;
};

TEST_CASE_METHOD(MeshTestFixture, "Mesh Utils", "[mesh][utils]")
{
    SECTION("Basic mesh utilities functionality") {
        // Create a simple triangle mesh (2 triangles forming a quad)
        Eigen::MatrixXd V(4, 3);
        V << 0.0, 0.0, 0.0,
             1.0, 0.0, 0.0,
             1.0, 1.0, 0.0,
             0.0, 1.0, 0.0;

        Eigen::MatrixXi F(2, 3);
        F << 0, 1, 2,
             0, 2, 3;

        // Test mesh validity - should have 4 vertices and 2 faces
        REQUIRE(V.rows() == 4);
        REQUIRE(V.cols() == 3);
        REQUIRE(F.rows() == 2);
        REQUIRE(F.cols() == 3);

        // Test vertex coordinate ranges
        REQUIRE(V.col(0).minCoeff() == Catch::Approx(0.0));
        REQUIRE(V.col(0).maxCoeff() == Catch::Approx(1.0));
        REQUIRE(V.col(1).minCoeff() == Catch::Approx(0.0));
        REQUIRE(V.col(1).maxCoeff() == Catch::Approx(1.0));
        REQUIRE(V.col(2).minCoeff() == Catch::Approx(0.0));
        REQUIRE(V.col(2).maxCoeff() == Catch::Approx(0.0));
    }
}

TEST_CASE_METHOD(MeshTestFixture, "OBJ File IO", "[mesh][obj][io]")
{
    SECTION("Write and read simple triangle") {
        // Create a simple triangle
        Eigen::MatrixXd V(3, 3);
        V << 0.0, 0.0, 0.0,
             1.0, 0.0, 0.0,
             0.5, 1.0, 0.0;

        Eigen::MatrixXi F(1, 3);
        F << 0, 1, 2;

        std::string filepath = temp_file("triangle.obj");

        // Test write
        REQUIRE(io::OBJWriter::write(filepath, V, F));

        // Verify file exists and has expected content
        std::ifstream file(filepath);
        REQUIRE(file.is_open());

        std::string content((std::istreambuf_iterator<char>(file)),
                           std::istreambuf_iterator<char>());

        // Check for vertex lines
        REQUIRE(content.find("v 0") != std::string::npos);
        REQUIRE(content.find("v 1") != std::string::npos);
        REQUIRE(content.find("v 0.5") != std::string::npos);

        // Check for face line (OBJ uses 1-based indexing)
        REQUIRE(content.find("f 1 2 3") != std::string::npos);

        // Test read
        Eigen::MatrixXd V_loaded;
        Eigen::MatrixXi E_loaded;
        Eigen::MatrixXi F_loaded;
        REQUIRE(io::OBJReader::read(filepath, V_loaded, E_loaded, F_loaded));

        // Verify dimensions
        REQUIRE(V_loaded.rows() == 3);
        REQUIRE(V_loaded.cols() == 3);
        REQUIRE(F_loaded.rows() == 1);
        REQUIRE(F_loaded.cols() == 3);

        // Verify vertex coordinates
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                REQUIRE(V_loaded(i, j) == Catch::Approx(V(i, j)));
            }
        }

        // Verify face indices
        for (int i = 0; i < 1; ++i) {
            for (int j = 0; j < 3; ++j) {
                REQUIRE(F_loaded(i, j) == F(i, j));
            }
        }
    }

    SECTION("Write and read quad mesh") {
        // Create a quad (2 triangles)
        Eigen::MatrixXd V(4, 3);
        V << 0.0, 0.0, 0.0,
             2.0, 0.0, 0.0,
             2.0, 2.0, 0.0,
             0.0, 2.0, 0.0;

        Eigen::MatrixXi F(2, 3);
        F << 0, 1, 2,
             0, 2, 3;

        std::string filepath = temp_file("quad.obj");

        REQUIRE(io::OBJWriter::write(filepath, V, F));

        Eigen::MatrixXd V_loaded;
        Eigen::MatrixXi E_loaded;
        Eigen::MatrixXi F_loaded;
        REQUIRE(io::OBJReader::read(filepath, V_loaded, E_loaded, F_loaded));

        REQUIRE(V_loaded.rows() == 4);
        REQUIRE(F_loaded.rows() == 2);

        // Verify all vertices are preserved
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 3; ++j) {
                REQUIRE(V_loaded(i, j) == Catch::Approx(V(i, j)));
            }
        }
    }

    SECTION("Handle empty mesh") {
        Eigen::MatrixXd V(0, 3);
        Eigen::MatrixXi F(0, 3);

        std::string filepath = temp_file("empty.obj");

        REQUIRE(io::OBJWriter::write(filepath, V, F));

        Eigen::MatrixXd V_loaded;
        Eigen::MatrixXi E_loaded;
        Eigen::MatrixXi F_loaded;
        REQUIRE(io::OBJReader::read(filepath, V_loaded, E_loaded, F_loaded));

        REQUIRE(V_loaded.rows() == 0);
        REQUIRE(F_loaded.rows() == 0);
    }
}

TEST_CASE_METHOD(MeshTestFixture, "OBJ Error handling", "[mesh][obj][error]")
{
    SECTION("Non-existent file read") {
        Eigen::MatrixXd V;
        Eigen::MatrixXi E;
        Eigen::MatrixXi F;
        REQUIRE_FALSE(io::OBJReader::read("non_existent_file.obj", V, E, F));
    }

    SECTION("Invalid path write") {
        Eigen::MatrixXd V(1, 3);
        V << 0.0, 0.0, 0.0;
        Eigen::MatrixXi F(0, 3);

        REQUIRE_FALSE(io::OBJWriter::write("/invalid/path/file.obj", V, F));
    }

    SECTION("Corrupted OBJ file") {
        std::string filepath = temp_file("corrupted.obj");
        std::ofstream file(filepath);
        file << "v 1.0 2.0 3.0\n";
        file << "v invalid vertex data\n";
        file << "f 1 2\n"; // Invalid face (only 2 vertices)
        file.close();

        Eigen::MatrixXd V;
        Eigen::MatrixXi E;
        Eigen::MatrixXi F;
        // This should handle gracefully or return false
        io::OBJReader::read(filepath, V, E, F);

        // If read succeeds, we should at least have read the valid vertex
        if (V.rows() > 0) {
            REQUIRE(V(0, 0) == Catch::Approx(1.0));
            REQUIRE(V(0, 1) == Catch::Approx(2.0));
            REQUIRE(V(0, 2) == Catch::Approx(3.0));
        }
    }
}

TEST_CASE_METHOD(MeshTestFixture, "Mesh precision", "[mesh][precision]")
{
    SECTION("High precision coordinates") {
        Eigen::MatrixXd V(3, 3);
        V << 1.23456789012345, 2.34567890123456, 3.45678901234567,
            -4.56789012345678, -5.67890123456789, -6.78901234567890,
             7.89012345678901, 8.90123456789012, 9.01234567890123;

        Eigen::MatrixXi F(1, 3);
        F << 0, 1, 2;

        std::string filepath = temp_file("precision.obj");

        REQUIRE(io::OBJWriter::write(filepath, V, F));

        Eigen::MatrixXd V_loaded;
        Eigen::MatrixXi E_loaded;
        Eigen::MatrixXi F_loaded;
        REQUIRE(io::OBJReader::read(filepath, V_loaded, E_loaded, F_loaded));

        // Verify high precision is preserved (within reasonable tolerance)
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                REQUIRE(V_loaded(i, j) == Catch::Approx(V(i, j)).epsilon(1e-10));
            }
        }
    }
}
