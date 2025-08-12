#include "polyfem.h"
#include <string>
#include <cstring>
#include <memory>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

// Simple JSON implementation to avoid external dependencies
namespace simple_json {
    std::string escape_string(const std::string& str) {
        std::string result;
        for (char c : str) {
            switch (c) {
                case '"': result += "\\\""; break;
                case '\\': result += "\\\\"; break;
                case '\n': result += "\\n"; break;
                case '\r': result += "\\r"; break;
                case '\t': result += "\\t"; break;
                default: result += c; break;
            }
        }
        return result;
    }
}

// Simple 3D vector class
struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Vec3 operator-(const Vec3& other) const { return Vec3(x - other.x, y - other.y, z - other.z); }
    Vec3 cross(const Vec3& other) const {
        return Vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
    }
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
};

// Simple mesh structure
struct SimpleMesh {
    std::vector<Vec3> vertices;
    std::vector<std::vector<int>> faces;
    
    bool load_obj(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string prefix;
            iss >> prefix;
            
            if (prefix == "v") {
                double x, y, z;
                if (iss >> x >> y >> z) {
                    vertices.emplace_back(x, y, z);
                }
            } else if (prefix == "f") {
                std::vector<int> face;
                std::string vertex_data;
                while (iss >> vertex_data) {
                    // Handle vertex/texture/normal format (v/vt/vn)
                    size_t slash_pos = vertex_data.find('/');
                    std::string vertex_index_str = vertex_data.substr(0, slash_pos);
                    int vertex_index = std::stoi(vertex_index_str) - 1; // OBJ is 1-indexed
                    face.push_back(vertex_index);
                }
                if (face.size() >= 3) {
                    faces.push_back(face);
                }
            }
        }
        
        return !vertices.empty() && !faces.empty();
    }
    
    Vec3 get_bbox_min() const {
        if (vertices.empty()) return Vec3();
        Vec3 min_pt = vertices[0];
        for (const auto& v : vertices) {
            min_pt.x = std::min(min_pt.x, v.x);
            min_pt.y = std::min(min_pt.y, v.y);
            min_pt.z = std::min(min_pt.z, v.z);
        }
        return min_pt;
    }
    
    Vec3 get_bbox_max() const {
        if (vertices.empty()) return Vec3();
        Vec3 max_pt = vertices[0];
        for (const auto& v : vertices) {
            max_pt.x = std::max(max_pt.x, v.x);
            max_pt.y = std::max(max_pt.y, v.y);
            max_pt.z = std::max(max_pt.z, v.z);
        }
        return max_pt;
    }
    
    double calculate_surface_area() const {
        double total_area = 0.0;
        for (const auto& face : faces) {
            if (face.size() >= 3) {
                Vec3 v0 = vertices[face[0]];
                Vec3 v1 = vertices[face[1]];
                Vec3 v2 = vertices[face[2]];
                Vec3 cross_product = (v1 - v0).cross(v2 - v0);
                total_area += 0.5 * cross_product.norm();
            }
        }
        return total_area;
    }
};

UNIFEX_TERM simulate(UnifexEnv* env, UnifexPayload* config, char* output_path) {
    // Input validation
    if (!config || !output_path) {
        return simulate_result_error(env, "Invalid parameters");
    }
    
    if (strlen(output_path) == 0) {
        return simulate_result_error(env, "Output path cannot be empty");
    }
    
    try {
        // Extract JSON configuration from payload
        std::string config_json;
        if (config->data && config->size > 0) {
            config_json = std::string(reinterpret_cast<char*>(config->data), config->size);
        } else {
            return simulate_result_error(env, "Empty configuration payload");
        }
        
        // Create output directory if it doesn't exist
        std::string output_dir(output_path);
        
        // For now, we'll create a simple simulation result
        // In a full implementation, this would:
        // 1. Parse the JSON configuration
        // 2. Load garment and avatar meshes
        // 3. Set up PolyFEM simulation parameters
        // 4. Run the cloth fitting simulation
        // 5. Save results to output directory
        
        // Create a simple result JSON
        std::ostringstream result_stream;
        result_stream << "{";
        result_stream << "\"status\":\"completed\",";
        result_stream << "\"output_path\":\"" << simple_json::escape_string(output_dir) << "\",";
        result_stream << "\"config_size\":" << config->size << ",";
        result_stream << "\"simulation_time\":0.0,";
        result_stream << "\"iterations\":0,";
        result_stream << "\"message\":\"Simulation framework ready - full PolyFEM integration pending\"";
        result_stream << "}";
        
        std::string result = result_stream.str();
        return simulate_result_ok(env, result.c_str());
        
    } catch (const std::exception& e) {
        std::string error_msg = "Simulation error: " + std::string(e.what());
        return simulate_result_error(env, error_msg.c_str());
    }
}

UNIFEX_TERM validate_garment_mesh(UnifexEnv* env, char* mesh_path) {
    // Input validation
    if (!mesh_path || strlen(mesh_path) == 0) {
        return validate_garment_mesh_result_error(env, "Invalid mesh path");
    }
    
    try {
        // Load mesh using simple OBJ reader
        SimpleMesh mesh;
        if (!mesh.load_obj(std::string(mesh_path))) {
            return validate_garment_mesh_result_error(env, "Failed to read mesh file");
        }
        
        // Basic validation checks for garment meshes
        bool is_valid = true;
        std::string validation_errors;
        
        // Check if mesh has vertices and faces
        if (mesh.vertices.empty()) {
            is_valid = false;
            validation_errors += "No vertices found. ";
        }
        
        if (mesh.faces.empty()) {
            is_valid = false;
            validation_errors += "No faces found. ";
        }
        
        // Check for reasonable mesh size (garments should have some complexity)
        if (mesh.vertices.size() < 10) {
            is_valid = false;
            validation_errors += "Too few vertices for a garment mesh. ";
        }
        
        if (mesh.faces.size() < 10) {
            is_valid = false;
            validation_errors += "Too few faces for a garment mesh. ";
        }
        
        // Check for valid face indices
        for (const auto& face : mesh.faces) {
            for (int vertex_idx : face) {
                if (vertex_idx < 0 || vertex_idx >= static_cast<int>(mesh.vertices.size())) {
                    is_valid = false;
                    validation_errors += "Invalid face index found. ";
                    break;
                }
            }
            if (!is_valid) break;
        }
        
        // Check for degenerate faces (faces with zero area)
        int degenerate_count = 0;
        for (const auto& face : mesh.faces) {
            if (face.size() >= 3) {
                Vec3 v0 = mesh.vertices[face[0]];
                Vec3 v1 = mesh.vertices[face[1]];
                Vec3 v2 = mesh.vertices[face[2]];
                double area = 0.5 * (v1 - v0).cross(v2 - v0).norm();
                if (area < 1e-12) {
                    degenerate_count++;
                }
            }
        }
        
        if (degenerate_count > static_cast<int>(mesh.faces.size()) * 0.1) { // More than 10% degenerate faces
            is_valid = false;
            validation_errors += "Too many degenerate faces. ";
        }
        
        // Check bounding box (garment should have reasonable dimensions)
        Vec3 bbox_min = mesh.get_bbox_min();
        Vec3 bbox_max = mesh.get_bbox_max();
        Vec3 bbox_size = bbox_max - bbox_min;
        
        double max_dimension = std::max({bbox_size.x, bbox_size.y, bbox_size.z});
        if (max_dimension < 1e-6) {
            is_valid = false;
            validation_errors += "Mesh has zero or near-zero dimensions. ";
        }
        
        if (!is_valid) {
            return validate_garment_mesh_result_error(env, validation_errors.c_str());
        }
        
        return validate_garment_mesh_result_ok(env, is_valid);
        
    } catch (const std::exception& e) {
        std::string error_msg = "Error validating garment mesh: " + std::string(e.what());
        return validate_garment_mesh_result_error(env, error_msg.c_str());
    }
}

UNIFEX_TERM validate_avatar_mesh(UnifexEnv* env, char* mesh_path) {
    // Input validation
    if (!mesh_path || strlen(mesh_path) == 0) {
        return validate_avatar_mesh_result_error(env, "Invalid mesh path");
    }
    
    try {
        // Load mesh using simple OBJ reader
        SimpleMesh mesh;
        if (!mesh.load_obj(std::string(mesh_path))) {
            return validate_avatar_mesh_result_error(env, "Failed to read mesh file");
        }
        
        // Basic validation checks for avatar meshes
        bool is_valid = true;
        std::string validation_errors;
        
        // Check if mesh has vertices and faces
        if (mesh.vertices.empty()) {
            is_valid = false;
            validation_errors += "No vertices found. ";
        }
        
        if (mesh.faces.empty()) {
            is_valid = false;
            validation_errors += "No faces found. ";
        }
        
        // Avatar meshes should be more complex than garment meshes
        if (mesh.vertices.size() < 100) {
            is_valid = false;
            validation_errors += "Too few vertices for an avatar mesh. ";
        }
        
        if (mesh.faces.size() < 100) {
            is_valid = false;
            validation_errors += "Too few faces for an avatar mesh. ";
        }
        
        // Check for valid face indices
        for (const auto& face : mesh.faces) {
            for (int vertex_idx : face) {
                if (vertex_idx < 0 || vertex_idx >= static_cast<int>(mesh.vertices.size())) {
                    is_valid = false;
                    validation_errors += "Invalid face index found. ";
                    break;
                }
            }
            if (!is_valid) break;
        }
        
        // Check bounding box and proportions (avatar should be humanoid-like)
        Vec3 bbox_min = mesh.get_bbox_min();
        Vec3 bbox_max = mesh.get_bbox_max();
        Vec3 bbox_size = bbox_max - bbox_min;
        
        double max_dimension = std::max({bbox_size.x, bbox_size.y, bbox_size.z});
        if (max_dimension < 1e-6) {
            is_valid = false;
            validation_errors += "Mesh has zero or near-zero dimensions. ";
        }
        
        // Check for reasonable aspect ratios (relaxed for various avatar types)
        double height_to_width = bbox_size.y / std::max(bbox_size.x, bbox_size.z);
        if (height_to_width < 0.2 || height_to_width > 20.0) {
            is_valid = false;
            validation_errors += "Extreme aspect ratio for avatar mesh. ";
        }
        
        // Check for degenerate faces
        int degenerate_count = 0;
        for (const auto& face : mesh.faces) {
            if (face.size() >= 3) {
                Vec3 v0 = mesh.vertices[face[0]];
                Vec3 v1 = mesh.vertices[face[1]];
                Vec3 v2 = mesh.vertices[face[2]];
                double area = 0.5 * (v1 - v0).cross(v2 - v0).norm();
                if (area < 1e-12) {
                    degenerate_count++;
                }
            }
        }
        
        if (degenerate_count > static_cast<int>(mesh.faces.size()) * 0.05) { // More than 5% degenerate faces
            is_valid = false;
            validation_errors += "Too many degenerate faces. ";
        }
        
        // Check mesh complexity (avatars should have reasonable detail)
        if (mesh.vertices.size() > 100000) {
            // Very high poly count - might cause performance issues
            // For now we'll allow it, but could add warnings
        }
        
        if (!is_valid) {
            return validate_avatar_mesh_result_error(env, validation_errors.c_str());
        }
        
        return validate_avatar_mesh_result_ok(env, is_valid);
        
    } catch (const std::exception& e) {
        std::string error_msg = "Error validating avatar mesh: " + std::string(e.what());
        return validate_avatar_mesh_result_error(env, error_msg.c_str());
    }
}

UNIFEX_TERM load_garment_info(UnifexEnv* env, char* garment_path) {
    // Input validation
    if (!garment_path || strlen(garment_path) == 0) {
        return load_garment_info_result_error(env, "Invalid garment path");
    }
    
    try {
        // Load garment mesh using simple OBJ reader
        SimpleMesh mesh;
        if (!mesh.load_obj(std::string(garment_path))) {
            return load_garment_info_result_error(env, "Failed to read garment mesh file");
        }
        
        // Calculate mesh statistics
        int vertex_count = static_cast<int>(mesh.vertices.size());
        int face_count = static_cast<int>(mesh.faces.size());
        
        // Calculate bounding box
        Vec3 bbox_min = mesh.get_bbox_min();
        Vec3 bbox_max = mesh.get_bbox_max();
        Vec3 bbox_size = bbox_max - bbox_min;
        
        // Create JSON response with mesh information using simple JSON
        std::ostringstream json_stream;
        json_stream << "{";
        json_stream << "\"vertex_count\":" << vertex_count << ",";
        json_stream << "\"face_count\":" << face_count << ",";
        json_stream << "\"bounding_box\":{";
        json_stream << "\"min\":[" << bbox_min.x << "," << bbox_min.y << "," << bbox_min.z << "],";
        json_stream << "\"max\":[" << bbox_max.x << "," << bbox_max.y << "," << bbox_max.z << "],";
        json_stream << "\"size\":[" << bbox_size.x << "," << bbox_size.y << "," << bbox_size.z << "]";
        json_stream << "},";
        json_stream << "\"mesh_type\":\"garment\",";
        json_stream << "\"file_path\":\"" << simple_json::escape_string(std::string(garment_path)) << "\"";
        json_stream << "}";
        
        std::string info_str = json_stream.str();
        
        // Return JSON string directly as a binary
        return load_garment_info_result_ok(env, info_str.c_str());
        
    } catch (const std::exception& e) {
        std::string error_msg = "Error loading garment info: " + std::string(e.what());
        return load_garment_info_result_error(env, error_msg.c_str());
    }
}

UNIFEX_TERM load_avatar_info(UnifexEnv* env, char* avatar_path) {
    // Input validation
    if (!avatar_path || strlen(avatar_path) == 0) {
        return load_avatar_info_result_error(env, "Invalid avatar path");
    }
    
    try {
        // Load avatar mesh using simple OBJ reader
        SimpleMesh mesh;
        if (!mesh.load_obj(std::string(avatar_path))) {
            return load_avatar_info_result_error(env, "Failed to read avatar mesh file");
        }
        
        // Calculate mesh statistics
        int vertex_count = static_cast<int>(mesh.vertices.size());
        int face_count = static_cast<int>(mesh.faces.size());
        
        // Calculate bounding box
        Vec3 bbox_min = mesh.get_bbox_min();
        Vec3 bbox_max = mesh.get_bbox_max();
        Vec3 bbox_size = bbox_max - bbox_min;
        
        // Calculate additional avatar-specific metrics
        double surface_area = mesh.calculate_surface_area();
        
        // Check if this looks like a valid avatar mesh (basic heuristics)
        bool is_likely_avatar = (vertex_count > 1000 && face_count > 1000 && 
                                bbox_size.y > bbox_size.x && bbox_size.y > bbox_size.z);
        
        // Create JSON response with mesh information using simple JSON
        std::ostringstream json_stream;
        json_stream << "{";
        json_stream << "\"vertex_count\":" << vertex_count << ",";
        json_stream << "\"face_count\":" << face_count << ",";
        json_stream << "\"surface_area\":" << surface_area << ",";
        json_stream << "\"bounding_box\":{";
        json_stream << "\"min\":[" << bbox_min.x << "," << bbox_min.y << "," << bbox_min.z << "],";
        json_stream << "\"max\":[" << bbox_max.x << "," << bbox_max.y << "," << bbox_max.z << "],";
        json_stream << "\"size\":[" << bbox_size.x << "," << bbox_size.y << "," << bbox_size.z << "]";
        json_stream << "},";
        json_stream << "\"mesh_type\":\"avatar\",";
        json_stream << "\"is_likely_avatar\":" << (is_likely_avatar ? "true" : "false") << ",";
        json_stream << "\"file_path\":\"" << simple_json::escape_string(std::string(avatar_path)) << "\"";
        json_stream << "}";
        
        std::string info_str = json_stream.str();
        
        // Return JSON string directly as a binary
        return load_avatar_info_result_ok(env, info_str.c_str());
        
    } catch (const std::exception& e) {
        std::string error_msg = "Error loading avatar info: " + std::string(e.what());
        return load_avatar_info_result_error(env, error_msg.c_str());
    }
}
