#pragma once

namespace polyfem {
namespace autogen {

    /// @brief Gradient of the python function
    /// def mix_product(t0, t1, t2):
    ///     return numpy.cross(t0, t1).dot(t2) / sympy.sqrt(t0.dot(t0) * t1.dot(t1) * t2.dot(t2))
    void angle_gradient(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[9]);

    /// @brief Hessian of the python function
    /// def mix_product(t0, t1, t2):
    ///     return numpy.cross(t0, t1).dot(t2) / sympy.sqrt(t0.dot(t0) * t1.dot(t1) * t2.dot(t2))
    // dA is (81×1) flattened in column-major order
    void angle_hessian(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[81]);

    /// @brief Gradient of the python function
    /// def triangle_area(t0, t1, t2):
    ///     n = numpy.cross(t1 - t0, t2 - t0)
    ///     return norm(n) / 2
    void triangle_area_gradient(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[9]);

    /// @brief Hessian of the python function
    /// def triangle_area(t0, t1, t2):
    ///     n = numpy.cross(t1 - t0, t2 - t0)
    ///     return norm(n) / 2
    // dA is (81×1) flattened in column-major order
    void triangle_area_hessian(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[81]);

    /// @brief Gradient of the python function
    /// def normalized_dot(a,b,c,d,e,f):
    ///     n1 = numpy.cross((b - a), (c - a))
    ///     n2 = numpy.cross((e - d), (f - d))
    ///     return n1.dot(n2) / sympy.sqrt(n1.dot(n1) * n2.dot(n2))
    void normal_dot_product_gradient(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double t3_x,
        double t3_y,
        double t3_z,
        double dA[12]);

    /// @brief Hessian of the python function
    /// def normalized_dot(a,b,c,d,e,f):
    ///     n1 = numpy.cross((b - a), (c - a))
    ///     n2 = numpy.cross((e - d), (f - d))
    ///     return n1.dot(n2) / sympy.sqrt(n1.dot(n1) * n2.dot(n2))
    // dA is (324×1) flattened in column-major order
    void normal_dot_product_hessian(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double t3_x,
        double t3_y,
        double t3_z,
        double dA[144]);

    void normal_triple_product_gradient(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double t3_x,
        double t3_y,
        double t3_z,
        double dA[12]);
    // dA is (144×1) flattened in column-major order
    void normal_triple_product_hessian(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double t3_x,
        double t3_y,
        double t3_z,
        double dA[144]);

    void point_edge_distance_gradient(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[9]);
    // dA is (81×1) flattened in column-major order
    void point_edge_distance_hessian(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[81]);


    void curve_cross_product_norm_gradient(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[9]);
    // dA is (81×1) flattened in column-major order
    void curve_cross_product_norm_hessian(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[81]);


    void curve_dot_product_norm_gradient(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[9]);
    // dA is (81×1) flattened in column-major order
    void curve_dot_product_norm_hessian(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[81]);
}
}