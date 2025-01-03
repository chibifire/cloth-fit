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


    void curve_dot_product_norm_gradient(
        double p1x,
        double p1y,
        double p1z,
        double p2x,
        double p2y,
        double p2z,
        double p3x,
        double p3y,
        double p3z,
        double grad[9]);
    // hess is (81×1) flattened in column-major order
    void curve_dot_product_norm_hessian(
        double p1x,
        double p1y,
        double p1z,
        double p2x,
        double p2y,
        double p2z,
        double p3x,
        double p3y,
        double p3z,
        double hess[81]);

    void curve_twist_gradient(double t0_x, double t0_y, double t0_z, double t1_x, double t1_y, double t1_z, double t2_x, double t2_y, double t2_z, double t3_x, double t3_y, double t3_z, double dA[12]);
    // dA is (144×1) flattened in column-major order
    void curve_twist_hessian(double t0_x, double t0_y, double t0_z, double t1_x, double t1_y, double t1_z, double t2_x, double t2_y, double t2_z, double t3_x, double t3_y, double t3_z, double dA[144]);

    void def_grad_gradient(
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
    void def_grad_hessian(
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

    void line_projection_uv_gradient(
        double t, double vx, double vy, double vz, double dA[4],
        double a0x, double a0y, double a0z,
        double a1x, double a1y, double a1z,
        double b0x, double b0y, double b0z,
        double b1x, double b1y, double b1z);
    // dA is (16×1) flattened in column-major order
    void line_projection_uv_hessian(
        double t, double vx, double vy, double vz, double dA[16],
        double a0x, double a0y, double a0z,
        double a1x, double a1y, double a1z,
        double b0x, double b0y, double b0z,
        double b1x, double b1y, double b1z);


    void curve_torsion_sin_gradient(
        double p0x,
        double p0y,
        double p0z,
        double p1x,
        double p1y,
        double p1z,
        double p2x,
        double p2y,
        double p2z,
        double p3x,
        double p3y,
        double p3z,
        double grad[12]);
    // hess is (144×1) flattened in column-major order
    void curve_torsion_sin_hessian(
        double p0x,
        double p0y,
        double p0z,
        double p1x,
        double p1y,
        double p1z,
        double p2x,
        double p2y,
        double p2z,
        double p3x,
        double p3y,
        double p3z,
        double hess[144]);

    void curve_torsion_cos_gradient(
        double p0x,
        double p0y,
        double p0z,
        double p1x,
        double p1y,
        double p1z,
        double p2x,
        double p2y,
        double p2z,
        double p3x,
        double p3y,
        double p3z,
        double grad[12]);
    // hess is (144×1) flattened in column-major order
    void curve_torsion_cos_hessian(
        double p0x,
        double p0y,
        double p0z,
        double p1x,
        double p1y,
        double p1z,
        double p2x,
        double p2y,
        double p2z,
        double p3x,
        double p3y,
        double p3z,
        double hess[144]);
}
}