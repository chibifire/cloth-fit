////////////////////////////////////////////////////////////////////////////////
#include <polyfem/io/GLTFReader.hpp>
#include <polyfem/io/GLTFWriter.hpp>
#include <polyfem/io/OBJReader.hpp>
#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/mesh/MeshUtils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>

#include <filesystem>
#include <iostream>
#include <memory>
////////////////////////////////////////////////////////////////////////////////

using namespace polyfem;
using namespace polyfem::io;

TEST_CASE("GLTF Round Trip - Simple Triangle", "[gltf][io]")
{
    // Create a simple triangle mesh
    std::vector<std::vector<double>> V = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.5, 1.0, 0.0}
    };
    
    std::vector<std::vector<int>> F = {
        {0, 1, 2}
    };
    
    std::vector<std::vector<int>> L; // Empty lines
    
    const std::string test_file = std::string(POLYFEM_TEST_DIR) + "/test_triangle.gltf";
    
    // Convert to Eigen matrices
    Eigen::MatrixXd V_eigen(V.size(), 3);
    for (size_t i = 0; i < V.size(); ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            V_eigen(i, j) = V[i][j];
        }
    }
    
    Eigen::MatrixXi F_eigen(F.size(), 3);
    for (size_t i = 0; i < F.size(); ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            F_eigen(i, j) = F[i][j];
        }
    }
    
    // Test writing
    REQUIRE(GLTFWriter::write(test_file, V_eigen, F_eigen));
    
    // Test reading back
    std::vector<std::vector<double>> V_read, TC_read, N_read;
    std::vector<std::vector<int>> F_read, FTC_read, FN_read, L_read;
    
    REQUIRE(GLTFReader::read(test_file, V_read, TC_read, N_read, F_read, FTC_read, FN_read, L_read));
    
    // Verify vertex data
    REQUIRE(V_read.size() == V.size());
    for (size_t i = 0; i < V.size(); ++i)
    {
        REQUIRE(V_read[i].size() == V[i].size());
        for (size_t j = 0; j < V[i].size(); ++j)
        {
            REQUIRE(V_read[i][j] == Catch::Approx(V[i][j]).epsilon(1e-6));
        }
    }
    
    // Verify face data
    REQUIRE(F_read.size() == F.size());
    for (size_t i = 0; i < F.size(); ++i)
    {
        REQUIRE(F_read[i].size() == F[i].size());
        for (size_t j = 0; j < F[i].size(); ++j)
        {
            REQUIRE(F_read[i][j] == F[i][j]);
        }
    }
    
    // Clean up
    std::filesystem::remove(test_file);
}

TEST_CASE("GLTF Eigen Matrix Interface", "[gltf][io][eigen]")
{
    // Create test data using Eigen matrices
    Eigen::MatrixXd V(4, 3);
    V << 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0,
         1.0, 1.0, 0.0,
         0.0, 1.0, 0.0;
    
    Eigen::MatrixXi F(2, 3);
    F << 0, 1, 2,
         0, 2, 3;
    
    const std::string test_file = std::string(POLYFEM_TEST_DIR) + "/test_quad.gltf";
    
    // Test writing using Eigen interface
    REQUIRE(GLTFWriter::write(test_file, V, F));
    
    // Test reading back using Eigen interface
    Eigen::MatrixXd V_read;
    Eigen::MatrixXi E_read, F_read;
    
    REQUIRE(GLTFReader::read(test_file, V_read, E_read, F_read));
    
    // Verify dimensions
    REQUIRE(V_read.rows() == V.rows());
    REQUIRE(V_read.cols() == V.cols());
    REQUIRE(F_read.rows() == F.rows());
    REQUIRE(F_read.cols() == F.cols());
    
    // Verify vertex data
    for (int i = 0; i < V.rows(); ++i)
    {
        for (int j = 0; j < V.cols(); ++j)
        {
            REQUIRE(V_read(i, j) == Catch::Approx(V(i, j)).epsilon(1e-6));
        }
    }
    
    // Verify face data
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int j = 0; j < F.cols(); ++j)
        {
            REQUIRE(F_read(i, j) == F(i, j));
        }
    }
    
    // Clean up
    std::filesystem::remove(test_file);
}

TEST_CASE("GLTF with Normals and Texture Coordinates", "[gltf][io][attributes]")
{
    // Create test data with normals and texture coordinates
    std::vector<std::vector<double>> V = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.5, 1.0, 0.0}
    };
    
    std::vector<std::vector<double>> N = {
        {0.0, 0.0, 1.0},
        {0.0, 0.0, 1.0},
        {0.0, 0.0, 1.0}
    };
    
    std::vector<std::vector<double>> TC = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.5, 1.0}
    };
    
    std::vector<std::vector<int>> F = {{0, 1, 2}};
    std::vector<std::vector<int>> FN = {{0, 1, 2}};
    std::vector<std::vector<int>> FTC = {{0, 1, 2}};
    std::vector<std::vector<int>> L;
    
    const std::string test_file = std::string(POLYFEM_TEST_DIR) + "/test_attributes.gltf";
    
    // Test writing with attributes
    REQUIRE(GLTFWriter::write(test_file, V, N, TC, F, FN, FTC, L));
    
    // Test reading back
    std::vector<std::vector<double>> V_read, TC_read, N_read;
    std::vector<std::vector<int>> F_read, FTC_read, FN_read, L_read;
    
    REQUIRE(GLTFReader::read(test_file, V_read, TC_read, N_read, F_read, FTC_read, FN_read, L_read));
    
    // Verify normals
    REQUIRE(N_read.size() == N.size());
    for (size_t i = 0; i < N.size(); ++i)
    {
        REQUIRE(N_read[i].size() == N[i].size());
        for (size_t j = 0; j < N[i].size(); ++j)
        {
            REQUIRE(N_read[i][j] == Catch::Approx(N[i][j]).epsilon(1e-6));
        }
    }
    
    // Verify texture coordinates
    REQUIRE(TC_read.size() == TC.size());
    for (size_t i = 0; i < TC.size(); ++i)
    {
        REQUIRE(TC_read[i].size() == TC[i].size());
        for (size_t j = 0; j < TC[i].size(); ++j)
        {
            REQUIRE(TC_read[i][j] == Catch::Approx(TC[i][j]).epsilon(1e-6));
        }
    }
    
    // Clean up
    std::filesystem::remove(test_file);
}

TEST_CASE("GLTF Binary Format", "[gltf][io][binary]")
{
    // Create test data
    Eigen::MatrixXd V(3, 3);
    V << 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0,
         0.5, 1.0, 0.0;
    
    Eigen::MatrixXi F(1, 3);
    F << 0, 1, 2;
    
    const std::string test_file = std::string(POLYFEM_TEST_DIR) + "/test_binary.glb";
    
    // Test writing binary GLTF
    REQUIRE(GLTFWriter::write(test_file, V, F));
    
    // Test reading back
    Eigen::MatrixXd V_read;
    Eigen::MatrixXi E_read, F_read;
    
    REQUIRE(GLTFReader::read(test_file, V_read, E_read, F_read));
    
    // Verify data integrity
    REQUIRE(V_read.rows() == V.rows());
    REQUIRE(F_read.rows() == F.rows());
    
    for (int i = 0; i < V.rows(); ++i)
    {
        for (int j = 0; j < V.cols(); ++j)
        {
            REQUIRE(V_read(i, j) == Catch::Approx(V(i, j)).epsilon(1e-6));
        }
    }
    
    // Clean up
    std::filesystem::remove(test_file);
}

TEST_CASE("GLTF Coordinate System Preservation", "[gltf][io][coordinates]")
{
    // Create a test mesh with clear directional information
    // This tests that GLTF coordinate system conventions are preserved
    std::vector<std::vector<double>> V = {
        {0.0, 0.0, 0.0},  // Origin
        {1.0, 0.0, 0.0},  // +X direction
        {0.0, 1.0, 0.0},  // +Y direction (up in GLTF)
        {0.0, 0.0, 1.0}   // +Z direction
    };
    
    std::vector<std::vector<int>> F = {
        {0, 1, 2}, // Triangle in XY plane
        {0, 2, 3}  // Triangle in YZ plane
    };
    
    std::vector<std::vector<int>> L;
    
    const std::string test_file = std::string(POLYFEM_TEST_DIR) + "/test_coordinates.gltf";
    
    // Write and read back
    REQUIRE(GLTFWriter::write(test_file, V, F, L));
    
    std::vector<std::vector<double>> V_read, TC_read, N_read;
    std::vector<std::vector<int>> F_read, FTC_read, FN_read, L_read;
    
    REQUIRE(GLTFReader::read(test_file, V_read, TC_read, N_read, F_read, FTC_read, FN_read, L_read));
    
    // Verify that coordinate system is preserved exactly (no rotation applied)
    REQUIRE(V_read.size() == V.size());
    for (size_t i = 0; i < V.size(); ++i)
    {
        REQUIRE(V_read[i].size() == 3);
        REQUIRE(V_read[i][0] == Catch::Approx(V[i][0]).epsilon(1e-12)); // X
        REQUIRE(V_read[i][1] == Catch::Approx(V[i][1]).epsilon(1e-12)); // Y
        REQUIRE(V_read[i][2] == Catch::Approx(V[i][2]).epsilon(1e-12)); // Z
    }
    
    // Verify that Y-up convention is maintained
    // The +Y vertex should still be at (0, 1, 0)
    REQUIRE(V_read[2][0] == Catch::Approx(0.0).epsilon(1e-12));
    REQUIRE(V_read[2][1] == Catch::Approx(1.0).epsilon(1e-12));
    REQUIRE(V_read[2][2] == Catch::Approx(0.0).epsilon(1e-12));
    
    // Clean up
    std::filesystem::remove(test_file);
}

TEST_CASE("GLTF Error Handling", "[gltf][io][error]")
{
    std::vector<std::vector<double>> V_read, TC_read, N_read;
    std::vector<std::vector<int>> F_read, FTC_read, FN_read, L_read;
    
    // Test reading non-existent file
    REQUIRE_FALSE(GLTFReader::read("/non/existent/file.gltf", V_read, TC_read, N_read, F_read, FTC_read, FN_read, L_read));
    
    // Test reading invalid file
    const std::string invalid_file = std::string(POLYFEM_TEST_DIR) + "/invalid.gltf";
    {
        std::ofstream file(invalid_file);
        file << "invalid json content";
    }
    
    REQUIRE_FALSE(GLTFReader::read(invalid_file, V_read, TC_read, N_read, F_read, FTC_read, FN_read, L_read));
    
    // Clean up
    std::filesystem::remove(invalid_file);
}

TEST_CASE("GLTF vs OBJ Comparison", "[gltf][obj][comparison]")
{
    // Test that existing OBJ files can be converted to GLTF and maintain data integrity
    const std::string obj_file = std::string(POLYFEM_TEST_DIR) + "/garment.obj";
    
    // Skip test if OBJ file doesn't exist
    if (!std::filesystem::exists(obj_file))
    {
        SKIP("Test OBJ file not found");
    }
    
    // Read OBJ file
    std::vector<std::vector<double>> V_obj, TC_obj, N_obj;
    std::vector<std::vector<int>> F_obj, FTC_obj, FN_obj, L_obj;
    
    REQUIRE(OBJReader::read(obj_file, V_obj, TC_obj, N_obj, F_obj, FTC_obj, FN_obj, L_obj));
    
    // Write as GLTF
    const std::string gltf_file = std::string(POLYFEM_TEST_DIR) + "/converted_garment.gltf";
    REQUIRE(GLTFWriter::write(gltf_file, V_obj, N_obj, TC_obj, F_obj, FN_obj, FTC_obj, L_obj));
    
    // Read back from GLTF
    std::vector<std::vector<double>> V_gltf, TC_gltf, N_gltf;
    std::vector<std::vector<int>> F_gltf, FTC_gltf, FN_gltf, L_gltf;
    
    REQUIRE(GLTFReader::read(gltf_file, V_gltf, TC_gltf, N_gltf, F_gltf, FTC_gltf, FN_gltf, L_gltf));
    
    // Compare vertex counts
    REQUIRE(V_gltf.size() == V_obj.size());
    REQUIRE(F_gltf.size() == F_obj.size());
    
    // Compare vertex data (allowing for small floating point differences)
    for (size_t i = 0; i < std::min(V_obj.size(), size_t(10)); ++i) // Test first 10 vertices
    {
        REQUIRE(V_gltf[i].size() == V_obj[i].size());
        for (size_t j = 0; j < V_obj[i].size(); ++j)
        {
            REQUIRE(V_gltf[i][j] == Catch::Approx(V_obj[i][j]).epsilon(1e-5));
        }
    }
    
    // Clean up
    std::filesystem::remove(gltf_file);
}

TEST_CASE("GLTF Large Mesh Performance", "[gltf][io][performance]")
{
    // Create a larger mesh to test performance and robustness
    const int grid_size = 10;
    std::vector<std::vector<double>> V;
    std::vector<std::vector<int>> F;
    
    // Generate grid vertices
    for (int i = 0; i <= grid_size; ++i)
    {
        for (int j = 0; j <= grid_size; ++j)
        {
            V.push_back({static_cast<double>(i), static_cast<double>(j), 0.0});
        }
    }
    
    // Generate grid faces
    for (int i = 0; i < grid_size; ++i)
    {
        for (int j = 0; j < grid_size; ++j)
        {
            int v00 = i * (grid_size + 1) + j;
            int v01 = i * (grid_size + 1) + (j + 1);
            int v10 = (i + 1) * (grid_size + 1) + j;
            int v11 = (i + 1) * (grid_size + 1) + (j + 1);
            
            // Two triangles per grid cell
            F.push_back({v00, v01, v10});
            F.push_back({v01, v11, v10});
        }
    }
    
    std::vector<std::vector<int>> L; // Empty lines
    
    const std::string test_file = std::string(POLYFEM_TEST_DIR) + "/test_large.gltf";
    
    // Test writing and reading large mesh
    REQUIRE(GLTFWriter::write(test_file, V, F, L));
    
    std::vector<std::vector<double>> V_read, TC_read, N_read;
    std::vector<std::vector<int>> F_read, FTC_read, FN_read, L_read;
    
    REQUIRE(GLTFReader::read(test_file, V_read, TC_read, N_read, F_read, FTC_read, FN_read, L_read));
    
    // Verify counts
    REQUIRE(V_read.size() == V.size());
    REQUIRE(F_read.size() == F.size());
    
    // Verify some sample vertices
    for (size_t i = 0; i < std::min(V.size(), size_t(5)); ++i)
    {
        REQUIRE(V_read[i].size() == V[i].size());
        for (size_t j = 0; j < V[i].size(); ++j)
        {
            REQUIRE(V_read[i][j] == Catch::Approx(V[i][j]).epsilon(1e-6));
        }
    }
    
    // Clean up
    std::filesystem::remove(test_file);
}
