#include "auto_derivatives.hpp"
#include <cmath>

namespace polyfem::autogen {

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
        double dA[9])
    {
        const auto t0 =
            std::pow(t0_x, 2) + std::pow(t0_y, 2) + std::pow(t0_z, 2);
        const auto t1 =
            std::pow(t1_x, 2) + std::pow(t1_y, 2) + std::pow(t1_z, 2);
        const auto t2 =
            std::pow(t2_x, 2) + std::pow(t2_y, 2) + std::pow(t2_z, 2);
        const auto t3 = std::pow(t0 * t1 * t2, -1.0 / 2.0);
        const auto t4 = t0_z * t1_y;
        const auto t5 = t0_y * t1_x;
        const auto t6 = t0_x * t1_z - t0_z * t1_x;
        const auto t7 =
            t2_x * (t0_y * t1_z - t4) - t2_y * t6 + t2_z * (t0_x * t1_y - t5);
        const auto t8 = t7 / t0;
        const auto t9 = t7 / t1;
        const auto t10 = t7 / t2;
        dA[0] = t3 * (-t0_x * t8 + t1_y * t2_z - t1_z * t2_y);
        dA[1] = -t3 * (t0_y * t8 + t1_x * t2_z - t1_z * t2_x);
        dA[2] = t3 * (-t0_z * t8 + t1_x * t2_y - t1_y * t2_x);
        dA[3] = -t3 * (t0_y * t2_z - t0_z * t2_y + t1_x * t9);
        dA[4] = t3 * (t0_x * t2_z - t0_z * t2_x - t1_y * t9);
        dA[5] = -t3 * (t0_x * t2_y - t0_y * t2_x + t1_z * t9);
        dA[6] = t3 * (t0_y * t1_z - t10 * t2_x - t4);
        dA[7] = -t3 * (t10 * t2_y + t6);
        dA[8] = t3 * (t0_x * t1_y - t10 * t2_z - t5);
    }

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
        double dA[81])
    {
        const auto t0 = t1_y * t2_z;
        const auto t1 = t1_z * t2_y;
        const auto t2 = t0 - t1;
        const auto t3 = std::pow(t0_x, 2);
        const auto t4 = std::pow(t0_y, 2);
        const auto t5 = std::pow(t0_z, 2);
        const auto t6 = t3 + t4 + t5;
        const auto t7 = 1.0 / t6;
        const auto t8 = t0_x * t1_z;
        const auto t9 = t0_z * t1_x;
        const auto t10 = t8 - t9;
        const auto t11 = t10 * t2_y;
        const auto t12 = -t11;
        const auto t13 = t0_y * t1_z;
        const auto t14 = -t0_z * t1_y + t13;
        const auto t15 = t14 * t2_x;
        const auto t16 = t0_x * t1_y;
        const auto t17 = t0_y * t1_x;
        const auto t18 = t16 - t17;
        const auto t19 = t18 * t2_z;
        const auto t20 = t15 + t19;
        const auto t21 = t12 + t20;
        const auto t22 = std::pow(t1_x, 2);
        const auto t23 = std::pow(t1_y, 2);
        const auto t24 = std::pow(t1_z, 2);
        const auto t25 = t22 + t23 + t24;
        const auto t26 = std::pow(t2_x, 2);
        const auto t27 = std::pow(t2_y, 2);
        const auto t28 = std::pow(t2_z, 2);
        const auto t29 = t26 + t27 + t28;
        const auto t30 = std::pow(t25 * t29 * t6, -1.0 / 2.0);
        const auto t31 = t30 * t7;
        const auto t32 = t1_z * t2_x;
        const auto t33 = t1_x * t2_z - t32;
        const auto t34 = t0_x * t7;
        const auto t35 = 3 * t21;
        const auto t36 = t31 * (t0_x * t33 - t0_y * t2 + t0_y * t34 * t35);
        const auto t37 = t1_x * t2_y;
        const auto t38 = -t1_y * t2_x + t37;
        const auto t39 =
            t31 * (3 * t0_x * t0_z * t21 * t7 - t0_x * t38 - t0_z * t2);
        const auto t40 = t0_y * t2_z;
        const auto t41 = t0_z * t2_y;
        const auto t42 = t40 - t41;
        const auto t43 = 1.0 / t25;
        const auto t44 = t2 * t43;
        const auto t45 = t1_x * t43;
        const auto t46 = t30 * (-t1_x * t44 + t21 * t34 * t45 + t34 * t42);
        const auto t47 = t0_x * t2_z;
        const auto t48 = -t0_z * t2_x + t47;
        const auto t49 = t21 * t7;
        const auto t50 = t43 * t49;
        const auto t51 = t30 * (t16 * t50 - t1_y * t44 + t2_z - t34 * t48);
        const auto t52 = t0_x * t2_y;
        const auto t53 = t0_y * t2_x;
        const auto t54 = t52 - t53;
        const auto t55 = t30 * (-t1_z * t44 - t2_y + t34 * t54 + t50 * t8);
        const auto t56 = 1.0 / t29;
        const auto t57 = t2 * t56;
        const auto t58 =
            t30 * (t0_x * t21 * t2_x * t56 * t7 - t14 * t34 - t2_x * t57);
        const auto t59 = t49 * t56;
        const auto t60 = t30 * (t10 * t34 - t1_z - t2_y * t57 + t52 * t59);
        const auto t61 = t30 * (-t18 * t34 + t1_y - t2_z * t57 + t47 * t59);
        const auto t62 = t11 - t19;
        const auto t63 = -t15 + t62;
        const auto t64 = t0_y * t7;
        const auto t65 = t31 * (-t0_y * t38 + t0_z * t33 + t0_z * t35 * t64);
        const auto t66 = t33 * t43;
        const auto t67 = t30 * (t17 * t50 + t1_x * t66 - t2_z + t42 * t64);
        const auto t68 = t1_y * t43;
        const auto t69 = t21 * t64;
        const auto t70 = t30 * (t1_y * t66 - t48 * t64 + t68 * t69);
        const auto t71 = t30 * (t13 * t50 + t1_z * t66 + t2_x + t54 * t64);
        const auto t72 = t33 * t56;
        const auto t73 = t30 * (-t14 * t64 + t1_z + t2_x * t72 + t53 * t59);
        const auto t74 = t2_y * t56;
        const auto t75 = t30 * (t10 * t64 + t2_y * t72 + t69 * t74);
        const auto t76 = t30 * (-t18 * t64 - t1_x + t2_z * t72 + t40 * t59);
        const auto t77 = t0_z * t7;
        const auto t78 = t38 * t43;
        const auto t79 = t30 * (-t1_x * t78 + t2_y + t42 * t77 + t50 * t9);
        const auto t80 = t30
            * (t0_z * t1_y * t21 * t43 * t7 - t1_y * t78 - t2_x - t48 * t77);
        const auto t81 = t1_z * t43;
        const auto t82 = t30 * (-t1_z * t78 + t21 * t77 * t81 + t54 * t77);
        const auto t83 = t38 * t56;
        const auto t84 = t30
            * (t0_z * t21 * t2_x * t56 * t7 - t14 * t77 - t1_y - t2_x * t83);
        const auto t85 = t30 * (t10 * t77 + t1_x - t2_y * t83 + t41 * t59);
        const auto t86 =
            t30 * (t0_z * t21 * t2_z * t56 * t7 - t18 * t77 - t2_z * t83);
        const auto t87 = t35 * t43;
        const auto t88 = t30 * t43;
        const auto t89 = t35 * t45;
        const auto t90 = t88 * (-t1_x * t48 + t1_y * t42 + t1_y * t89);
        const auto t91 = t88 * (t1_x * t54 + t1_z * t42 + t1_z * t89);
        const auto t92 = t42 * t56;
        const auto t93 =
            t30 * (-t14 * t45 + t21 * t2_x * t45 * t56 + t2_x * t92);
        const auto t94 = t21 * t43 * t56;
        const auto t95 = t30 * (t0_z + t10 * t45 + t2_y * t92 + t37 * t94);
        const auto t96 = t30
            * (-t0_y - t18 * t45 + t1_x * t21 * t2_z * t43 * t56
               + t2_z * t42 * t56);
        const auto t97 = t88 * (t1_y * t54 + t1_z * t35 * t68 - t1_z * t48);
        const auto t98 = t48 * t56;
        const auto t99 = t30
            * (-t0_z - t14 * t68 + t1_y * t21 * t2_x * t43 * t56 - t2_x * t98);
        const auto t100 = t30 * (t10 * t68 + t21 * t68 * t74 - t2_y * t98);
        const auto t101 = t30 * (t0 * t94 + t0_x - t18 * t68 - t2_z * t98);
        const auto t102 = t54 * t56;
        const auto t103 = t30 * (t0_y + t102 * t2_x - t14 * t81 + t32 * t94);
        const auto t104 = t30 * (-t0_x + t1 * t94 + t10 * t81 + t102 * t2_y);
        const auto t105 =
            t30 * (t102 * t2_z - t18 * t81 + t21 * t2_z * t56 * t81);
        const auto t106 = t35 * t56;
        const auto t107 = t30 * t56;
        const auto t108 = t107 * (t10 * t2_x + t106 * t2_x * t2_y - t14 * t2_y);
        const auto t109 =
            t107 * (-t14 * t2_z - t18 * t2_x + 3 * t21 * t2_x * t2_z * t56);
        const auto t110 = t107 * (t10 * t2_z + t106 * t2_y * t2_z - t18 * t2_y);
        dA[0] = t31 * (-2 * t0_x * t2 + 3 * t21 * t3 * t7 - t21);
        dA[1] = t36;
        dA[2] = t39;
        dA[3] = t46;
        dA[4] = t51;
        dA[5] = t55;
        dA[6] = t58;
        dA[7] = t60;
        dA[8] = t61;
        dA[9] = t36;
        dA[10] = t31 * (2 * t0_y * t33 + t35 * t4 * t7 + t63);
        dA[11] = t65;
        dA[12] = t67;
        dA[13] = t70;
        dA[14] = t71;
        dA[15] = t73;
        dA[16] = t75;
        dA[17] = t76;
        dA[18] = t39;
        dA[19] = t65;
        dA[20] = t31 * (-2 * t0_z * t38 + 3 * t21 * t5 * t7 - t21);
        dA[21] = t79;
        dA[22] = t80;
        dA[23] = t82;
        dA[24] = t84;
        dA[25] = t85;
        dA[26] = t86;
        dA[27] = t46;
        dA[28] = t67;
        dA[29] = t79;
        dA[30] = t88 * (2 * t1_x * t42 + t22 * t87 + t63);
        dA[31] = t90;
        dA[32] = t91;
        dA[33] = t93;
        dA[34] = t95;
        dA[35] = t96;
        dA[36] = t51;
        dA[37] = t70;
        dA[38] = t80;
        dA[39] = t90;
        dA[40] = t88 * (-2 * t1_y * t48 + 3 * t21 * t23 * t43 - t21);
        dA[41] = t97;
        dA[42] = t99;
        dA[43] = t100;
        dA[44] = t101;
        dA[45] = t55;
        dA[46] = t71;
        dA[47] = t82;
        dA[48] = t91;
        dA[49] = t97;
        dA[50] = t88 * (2 * t1_z * t54 + t24 * t87 + t63);
        dA[51] = t103;
        dA[52] = t104;
        dA[53] = t105;
        dA[54] = t58;
        dA[55] = t73;
        dA[56] = t84;
        dA[57] = t93;
        dA[58] = t99;
        dA[59] = t103;
        dA[60] = t107 * (t106 * t26 - 3 * t15 + t62);
        dA[61] = t108;
        dA[62] = t109;
        dA[63] = t60;
        dA[64] = t75;
        dA[65] = t85;
        dA[66] = t95;
        dA[67] = t100;
        dA[68] = t104;
        dA[69] = t108;
        dA[70] = t107 * (3 * t10 * t2_y - t20 + 3 * t21 * t27 * t56);
        dA[71] = t110;
        dA[72] = t61;
        dA[73] = t76;
        dA[74] = t86;
        dA[75] = t96;
        dA[76] = t101;
        dA[77] = t105;
        dA[78] = t109;
        dA[79] = t110;
        dA[80] = t107 * (-t12 - t15 - 3 * t19 + 3 * t21 * t28 * t56);
    }

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
        double dA[9])
    {
        const auto t0 = -t2_y;
        const auto t1 = t0 + t1_y;
        const auto t2 = t0_x - t1_x;
        const auto t3 = t0 + t0_y;
        const auto t4 = -t2_x;
        const auto t5 = t0_x + t4;
        const auto t6 = t0_y - t1_y;
        const auto t7 = t2 * t3 - t5 * t6;
        const auto t8 = -t2_z;
        const auto t9 = t1_z + t8;
        const auto t10 = t0_z + t8;
        const auto t11 = t0_z - t1_z;
        const auto t12 = t10 * t2 - t11 * t5;
        const auto t13 = t10 * t6 - t11 * t3;
        const auto t14 = (1.0 / 2.0)
            / std::sqrt(std::pow(t12, 2) + std::pow(t13, 2) + std::pow(t7, 2));
        const auto t15 = t1_x + t4;
        const auto t16 = t10 * t2 - t11 * t5;
        dA[0] = t14 * (t1 * t7 + t12 * t9);
        dA[1] = -t14 * (-t13 * t9 + t15 * t7);
        dA[2] = -t14 * (t1 * t13 + t15 * t16);
        dA[3] = -t14 * (t10 * t16 + t3 * t7);
        dA[4] = t14 * (-t10 * t13 + t5 * t7);
        dA[5] = t14 * (t12 * t5 + t13 * t3);
        dA[6] = t14 * (t11 * t12 + t6 * t7);
        dA[7] = -t14 * (-t11 * t13 + t2 * t7);
        dA[8] = -t14 * (t13 * t6 + t16 * t2);
    }

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
        double dA[81])
    {
        const auto t0 = -t2_y;
        const auto t1 = t0 + t1_y;
        const auto t2 = std::pow(t1, 2);
        const auto t3 = -t2_z;
        const auto t4 = t1_z + t3;
        const auto t5 = std::pow(t4, 2);
        const auto t6 = t0_x - t1_x;
        const auto t7 = t0 + t0_y;
        const auto t8 = t6 * t7;
        const auto t9 = -t2_x;
        const auto t10 = t0_x + t9;
        const auto t11 = t0_y - t1_y;
        const auto t12 = t10 * t11;
        const auto t13 = -t12;
        const auto t14 = t13 + t8;
        const auto t15 = t1 * t14;
        const auto t16 = t0_z + t3;
        const auto t17 = t16 * t6;
        const auto t18 = t0_z - t1_z;
        const auto t19 = t10 * t18;
        const auto t20 = -t19;
        const auto t21 = t17 + t20;
        const auto t22 = t15 + t21 * t4;
        const auto t23 = t11 * t16;
        const auto t24 = t18 * t7;
        const auto t25 = -t24;
        const auto t26 = t23 + t25;
        const auto t27 = std::pow(t14, 2) + std::pow(t26, 2);
        const auto t28 = std::pow(t21, 2) + t27;
        const auto t29 = 1.0 / t28;
        const auto t30 = -t6;
        const auto t31 = -t16;
        const auto t32 = -t10;
        const auto t33 = -t18;
        const auto t34 = t30 * t31 - t32 * t33;
        const auto t35 = t15 + t34 * t4;
        const auto t36 = t29 * t35;
        const auto t37 = (1.0 / 2.0) / std::sqrt(t28);
        const auto t38 = t1_x + t9;
        const auto t39 = t1 * t38;
        const auto t40 = t14 * t38 - t26 * t4;
        const auto t41 = -t40;
        const auto t42 = t27 + std::pow(t34, 2);
        const auto t43 = 1.0 / t42;
        const auto t44 = -t7;
        const auto t45 = -t11;
        const auto t46 = t30 * t44 - t32 * t45;
        const auto t47 = -t34;
        const auto t48 = t1 * t46 - t4 * t47;
        const auto t49 = t43 * t48;
        const auto t50 = (1.0 / 2.0) / std::sqrt(t42);
        const auto t51 = t38 * t4;
        const auto t52 = t1 * t26 + t34 * t38;
        const auto t53 = t14 * t7 + t16 * t34;
        const auto t54 = t1 * t7;
        const auto t55 = t16 * t4;
        const auto t56 = t54 + t55;
        const auto t57 = t10 * t14 - t16 * t26;
        const auto t58 = t1 * t10 + t14;
        const auto t59 = t26 * t7;
        const auto t60 = t10 * t34 + t59;
        const auto t61 = t10 * t4 + t21;
        const auto t62 = t11 * t14;
        const auto t63 = t18 * t21 + t62;
        const auto t64 = t1 * t11;
        const auto t65 = t18 * t4;
        const auto t66 = t64 + t65;
        const auto t67 = t14 * t6 - t18 * t26;
        const auto t68 = -t67;
        const auto t69 = t1 * t6 + t14;
        const auto t70 = t11 * t26 + t34 * t6;
        const auto t71 = t21 + t4 * t6;
        const auto t72 = t31 * t45 - t33 * t44;
        const auto t73 = -t38 * t46 + t4 * t72;
        const auto t74 = t43 * t73;
        const auto t75 = std::pow(t38, 2);
        const auto t76 = t1 * t4;
        const auto t77 = t12 - t8;
        const auto t78 = t29 * t41;
        const auto t79 = t10 * t38;
        const auto t80 = t55 + t79;
        const auto t81 = t26 + t4 * t7;
        const auto t82 = t18 * t34 + t62;
        const auto t83 = t11 * t38 + t77;
        const auto t84 = t38 * t6;
        const auto t85 = t65 + t84;
        const auto t86 = t11 * t4 + t26;
        const auto t87 = t43 * (-t1 * t72 + t38 * t47);
        const auto t88 = -t17 + t19;
        const auto t89 = t16 * t38 + t88;
        const auto t90 = t10 * t21 + t59;
        const auto t91 = t54 + t79;
        const auto t92 = t18 * t38;
        const auto t93 = t1 * t18;
        const auto t94 = -t23 + t24;
        const auto t95 = t37 * (-t29 * t52 * t70 + t64 + t84);
        const auto t96 = t43 * (t16 * t47 + t44 * t46);
        const auto t97 = std::pow(t7, 2);
        const auto t98 = std::pow(t16, 2);
        const auto t99 = t10 * t7;
        const auto t100 = t10 * t16;
        const auto t101 = t11 * t7;
        const auto t102 = t16 * t18;
        const auto t103 = t101 + t102;
        const auto t104 = 2 * t17 + t20;
        const auto t105 = t10 * t46 + t31 * t72;
        const auto t106 = t105 * t43;
        const auto t107 = std::pow(t10, 2);
        const auto t108 = t16 * t7;
        const auto t109 = -2 * t10 * t11 + t8;
        const auto t110 = t10 * t6;
        const auto t111 = t102 + t110;
        const auto t112 = t32 * t47 + t7 * t72;
        const auto t113 = t112 * t43;
        const auto t114 = -2 * t10 * t18 + t17;
        const auto t115 = -2 * t18 * t7 + t23;
        const auto t116 = t101 + t110;
        const auto t117 = t29 * t82;
        const auto t118 = t11 * t46 + t33 * t47;
        const auto t119 = t118 * t43;
        const auto t120 = std::pow(t11, 2);
        const auto t121 = std::pow(t18, 2);
        const auto t122 = t11 * t6;
        const auto t123 = t18 * t6;
        const auto t124 = t18 * t72 + t30 * t46;
        const auto t125 = t124 * t43;
        const auto t126 = std::pow(t6, 2);
        const auto t127 = t11 * t18;
        const auto t128 = t43 * (t45 * t72 + t47 * t6);
        dA[0] = t37 * (t2 - t22 * t36 + t5);
        dA[1] = -t50 * (t39 + t41 * t49);
        dA[2] = t50 * (t43 * t48 * t52 - t51);
        dA[3] = t37 * (t29 * t35 * t53 - t56);
        dA[4] = t50 * (-t49 * t57 + t58);
        dA[5] = t50 * (-t49 * t60 + t61);
        dA[6] = t37 * (-t36 * t63 + t66);
        dA[7] = t50 * (-t49 * t68 - t69);
        dA[8] = t50 * (t43 * t48 * t70 - t71);
        dA[9] = -t50 * (t35 * t74 + t39);
        dA[10] = t50 * (-t41 * t74 + t5 + t75);
        dA[11] = t50 * (t43 * t52 * t73 - t76);
        dA[12] = t50 * (t38 * t7 + t53 * t74 + t77);
        dA[13] = -t37 * (t57 * t78 + t80);
        dA[14] = t50 * (-t60 * t74 + t81);
        dA[15] = t50 * (-t74 * t82 - t83);
        dA[16] = t37 * (t67 * t78 + t85);
        dA[17] = t50 * (t43 * t70 * t73 - t86);
        dA[18] = -t50 * (t35 * t87 + t51);
        dA[19] = -t50 * (t41 * t87 + t76);
        dA[20] = t37 * (t2 - t29 * std::pow(t52, 2) + t75);
        dA[21] = t50 * (t53 * t87 + t89);
        dA[22] = t50 * (t1 * t16 - t26 - t57 * t87);
        dA[23] = t37 * (t29 * t52 * t90 - t91);
        dA[24] = t50 * (-t82 * t87 - t88 - t92);
        dA[25] = t50 * (-t68 * t87 - t93 - t94);
        dA[26] = t95;
        dA[27] = -t50 * (t35 * t96 + t56);
        dA[28] = t50 * (-t14 + t38 * t7 - t41 * t96);
        dA[29] = t50 * (t52 * t96 + t89);
        dA[30] = t50 * (t53 * t96 + t97 + t98);
        dA[31] = -t50 * (t57 * t96 + t99);
        dA[32] = -t50 * (t100 + t60 * t96);
        dA[33] = -t50 * (t103 + t82 * t96);
        dA[34] = t50 * (-t12 + 2 * t6 * t7 - t68 * t96);
        dA[35] = t50 * (t104 + t70 * t96);
        dA[36] = t50 * (-t106 * t35 + t58);
        dA[37] = t37 * (t29 * t40 * t57 - t80);
        dA[38] = t50 * (t1 * t16 + t106 * t52 + t94);
        dA[39] = t50 * (t105 * t43 * t53 - t99);
        dA[40] = t37 * (t107 - t29 * std::pow(t57, 2) + t98);
        dA[41] = -t50 * (t106 * t60 + t108);
        dA[42] = t50 * (-t106 * t82 - t109);
        dA[43] = t37 * (-t111 + t29 * t57 * t67);
        dA[44] = t50 * (t106 * t70 + 2 * t23 + t25);
        dA[45] = t50 * (-t113 * t35 + t61);
        dA[46] = t50 * (-t113 * t41 + t81);
        dA[47] = t50 * (t112 * t43 * t52 - t91);
        dA[48] = t50 * (-t100 + t112 * t43 * t53);
        dA[49] = -t50 * (t108 + t113 * t57);
        dA[50] = t50 * (t107 - t113 * t60 + t97);
        dA[51] = t50 * (-t113 * t82 - t114);
        dA[52] = t50 * (-t113 * t68 - t115);
        dA[53] = t50 * (t112 * t43 * t70 - t116);
        dA[54] = t37 * (-t117 * t22 + t66);
        dA[55] = t50 * (-t119 * t41 - t83);
        dA[56] = t50 * (t119 * t52 + t21 - t92);
        dA[57] = t37 * (-t103 + t29 * t53 * t82);
        dA[58] = t50 * (-t109 - t119 * t57);
        dA[59] = t50 * (-t114 - t119 * t60);
        dA[60] = t37 * (-t117 * t63 + t120 + t121);
        dA[61] = -t50 * (t119 * t68 + t122);
        dA[62] = t50 * (t118 * t43 * t70 - t123);
        dA[63] = t50 * (-t125 * t35 - t69);
        dA[64] = t50 * (-t125 * t41 + t85);
        dA[65] = t50 * (t125 * t52 + t26 - t93);
        dA[66] = t50 * (t125 * t53 + t13 + 2 * t8);
        dA[67] = -t50 * (t111 + t125 * t57);
        dA[68] = t50 * (-t115 - t125 * t60);
        dA[69] = -t50 * (t122 + t125 * t82);
        dA[70] = t50 * (t121 - t125 * t68 + t126);
        dA[71] = t50 * (t124 * t43 * t70 - t127);
        dA[72] = t50 * (-t128 * t35 - t71);
        dA[73] = t50 * (-t128 * t41 - t86);
        dA[74] = t95;
        dA[75] = t50 * (t104 + t128 * t53);
        dA[76] = t50 * (2 * t11 * t16 - t128 * t57 - t24);
        dA[77] = t37 * (-t116 + t29 * t70 * t90);
        dA[78] = -t50 * (t123 + t128 * t82);
        dA[79] = -t50 * (t127 + t128 * t68);
        dA[80] = t37 * (t120 + t126 - t29 * std::pow(t70, 2));
    }


    void normalized_dot_gradient(
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
        double dA[12])
    {
        const auto t0 = t0_x - t1_x;
        const auto t1 = -t0;
        const auto t2 = -t3_y;
        const auto t3 = t0_y + t2;
        const auto t4 = -t3;
        const auto t5 = -t3_x;
        const auto t6 = t0_x + t5;
        const auto t7 = -t6;
        const auto t8 = t0_y - t1_y;
        const auto t9 = -t8;
        const auto t10 = t1 * t4 - t7 * t9;
        const auto t11 = -t3_z;
        const auto t12 = t0_z + t11;
        const auto t13 = t0_z - t1_z;
        const auto t14 = t0 * t12 - t13 * t6;
        const auto t15 = std::pow(t14, 2);
        const auto t16 = -t12;
        const auto t17 = -t13;
        const auto t18 = t16 * t9 - t17 * t4;
        const auto t19 = std::pow(t10, 2) + t15 + std::pow(t18, 2);
        const auto t20 = -t2_z;
        const auto t21 = t0_z + t20;
        const auto t22 = -t21;
        const auto t23 = -t2_x;
        const auto t24 = t0_x + t23;
        const auto t25 = -t24;
        const auto t26 = t1 * t22 - t17 * t25;
        const auto t27 = -t2_y;
        const auto t28 = t0_y + t27;
        const auto t29 = t0 * t28 - t24 * t8;
        const auto t30 = -t13 * t28 + t21 * t8;
        const auto t31 = std::pow(t29, 2) + std::pow(t30, 2);
        const auto t32 = std::pow(t26, 2) + t31;
        const auto t33 = std::pow(t19 * t32, -1.0 / 2.0);
        const auto t34 = t1_y + t27;
        const auto t35 = t1_y + t2;
        const auto t36 = t1_z + t20;
        const auto t37 = t11 + t1_z;
        const auto t38 = 1.0 / t19;
        const auto t39 = 1.0 / t32;
        const auto t40 = t10 * t29 + t14 * t26 + t18 * t30;
        const auto t41 = -t10;
        const auto t42 = t1 * t16 - t17 * t7;
        const auto t43 = -t28;
        const auto t44 = t1 * t43 - t25 * t9;
        const auto t45 = -t26;
        const auto t46 = -t17 * t43 + t22 * t9;
        const auto t47 = std::pow(t44, 2) + std::pow(t45, 2) + std::pow(t46, 2);
        const auto t48 = -t18;
        const auto t49 = std::pow(t41, 2) + std::pow(t42, 2) + std::pow(t48, 2);
        const auto t50 = t1_x + t23;
        const auto t51 = t1_x + t5;
        const auto t52 = t39 * t40;
        const auto t53 = t38 * t52;
        const auto t54 = t0 * t21 - t13 * t24;
        const auto t55 = t31 + std::pow(t54, 2);
        const auto t56 = t12 * t8 - t13 * t3;
        const auto t57 =
            t15 + std::pow(t56, 2) + std::pow(t0 * t3 - t6 * t8, 2);
        const auto t58 = std::pow(t55 * t57, -1.0 / 2.0);
        const auto t59 = t29 * t8;
        const auto t60 = 1.0 / t55;
        const auto t61 = t10 * t8 + t13 * t14;
        const auto t62 = t0 * t29;
        const auto t63 = t13 * t30;
        const auto t64 = t0 * t10 - t13 * t18;
        const auto t65 = t0 * t14;
        const auto t66 = t0 * t26 + t30 * t8;
        dA[0] = t33
            * (-t10 * t34 - t14 * t36 - t26 * t37 - t29 * t35
               + t38 * t39 * t40
                   * (t47 * (-t35 * t41 + t37 * t42)
                      + t49 * (t34 * t44 - t36 * t45)));
        dA[1] = t33
            * (t10 * t50 - t18 * t36 + t29 * t51 - t30 * t37
               + t53
                   * (t47 * (-t37 * t48 + t41 * t51)
                      + t49 * (t36 * t46 - t44 * t50)));
        dA[2] = t33
            * (t14 * t50 + t18 * t34 + t26 * t51 + t30 * t35
               + t53
                   * (t47 * (t35 * t48 - t42 * t51)
                      + t49 * (-t34 * t46 + t45 * t50)));
        dA[3] = t33
            * (t10 * t28 + t12 * t26 + t14 * t21 + t29 * t3
               + t53
                   * (t47 * (t16 * t42 + t3 * t41)
                      + t49 * (t21 * t45 + t43 * t44)));
        dA[4] = t33
            * (-t10 * t24 + t12 * t30 + t18 * t21 - t29 * t6
               + t53
                   * (t47 * (t12 * t48 + t41 * t7)
                      + t49 * (t22 * t46 + t24 * t44)));
        dA[5] = t33
            * (-t14 * t24 - t18 * t28 - t26 * t6 - t3 * t30
               + t38 * t39 * t40
                   * (t47 * (t4 * t48 + t42 * t6)
                      + t49 * (t25 * t45 + t28 * t46)));
        dA[6] = t58 * (t40 * t60 * (t13 * t54 + t59) - t61);
        dA[7] = t33 * (t52 * (-t62 + t63) + t64);
        dA[8] = t58 * (-t40 * t60 * t66 + t56 * t8 + t65);
        dA[9] = t33 * (-t13 * t26 + t38 * t40 * t61 - t59);
        dA[10] = t58 * (-t40 * t64 / t57 + t62 - t63);
        dA[11] = t33 * (-t38 * t40 * (t18 * t8 + t65) + t66);
    }

    // dA is (144×1) flattened in column-major order
    void normalized_dot_hessian(
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
        double dA[144])
    {
        const auto t0 = t0_x - t1_x;
        const auto t1 = -t2_z;
        const auto t2 = t0_z + t1;
        const auto t3 = t0 * t2;
        const auto t4 = -t2_x;
        const auto t5 = t0_x + t4;
        const auto t6 = t0_z - t1_z;
        const auto t7 = t5 * t6;
        const auto t8 = -t7;
        const auto t9 = t3 + t8;
        const auto t10 = -t2_y;
        const auto t11 = t0_y + t10;
        const auto t12 = t0 * t11;
        const auto t13 = t0_y - t1_y;
        const auto t14 = t13 * t5;
        const auto t15 = -t14;
        const auto t16 = t12 + t15;
        const auto t17 = t13 * t2;
        const auto t18 = t11 * t6;
        const auto t19 = -t18;
        const auto t20 = t17 + t19;
        const auto t21 = std::pow(t16, 2) + std::pow(t20, 2);
        const auto t22 = t21 + std::pow(t9, 2);
        const auto t23 = -t3_y;
        const auto t24 = t0_y + t23;
        const auto t25 = t0 * t24;
        const auto t26 = -t3_x;
        const auto t27 = t0_x + t26;
        const auto t28 = t13 * t27;
        const auto t29 = -t28;
        const auto t30 = t25 + t29;
        const auto t31 = -t3_z;
        const auto t32 = t0_z + t31;
        const auto t33 = t0 * t32;
        const auto t34 = t27 * t6;
        const auto t35 = -t34;
        const auto t36 = t33 + t35;
        const auto t37 = std::pow(t36, 2);
        const auto t38 = t13 * t32;
        const auto t39 = t24 * t6;
        const auto t40 = -t39;
        const auto t41 = t38 + t40;
        const auto t42 = std::pow(t30, 2) + t37 + std::pow(t41, 2);
        const auto t43 = std::pow(t22 * t42, -1.0 / 2.0);
        const auto t44 = t10 + t1_y;
        const auto t45 = t1_y + t23;
        const auto t46 = 2 * t44 * t45;
        const auto t47 = t1 + t1_z;
        const auto t48 = t1_z + t31;
        const auto t49 = 2 * t47 * t48;
        const auto t50 = 1.0 / t22;
        const auto t51 = 1.0 / t42;
        const auto t52 = -t0;
        const auto t53 = -t24;
        const auto t54 = t52 * t53;
        const auto t55 = -t27;
        const auto t56 = -t13;
        const auto t57 = t55 * t56;
        const auto t58 = t54 - t57;
        const auto t59 = -t2;
        const auto t60 = t52 * t59;
        const auto t61 = -t5;
        const auto t62 = -t6;
        const auto t63 = t61 * t62;
        const auto t64 = t60 - t63;
        const auto t65 = t16 * t45 + t36 * t47 + t44 * t58 + t48 * t64;
        const auto t66 = t16 * t44 + t47 * t64;
        const auto t67 = -t32;
        const auto t68 = t56 * t67;
        const auto t69 = -t53 * t62 + t68;
        const auto t70 = t37 + std::pow(t58, 2) + std::pow(t69, 2);
        const auto t71 = t36 * t48 + t45 * t58;
        const auto t72 = t21 + std::pow(t64, 2);
        const auto t73 = t66 * t70 + t71 * t72;
        const auto t74 = std::pow(t42, -2);
        const auto t75 = t16 * t58 + t20 * t69 + t36 * t64;
        const auto t76 = std::pow(t22, -2);
        const auto t77 = t75 * t76;
        const auto t78 = t74 * t77;
        const auto t79 = std::pow(t44, 2);
        const auto t80 = std::pow(t47, 2);
        const auto t81 = std::pow(t45, 2);
        const auto t82 = std::pow(t48, 2);
        const auto t83 = 2 * t66;
        const auto t84 = t51 * t77;
        const auto t85 = t83 * t84;
        const auto t86 = 2 * t71;
        const auto t87 = t50 * t74 * t75;
        const auto t88 = t86 * t87;
        const auto t89 = t1_x + t4;
        const auto t90 = -t16 * t89 + t20 * t47;
        const auto t91 = t1_x + t26;
        const auto t92 = t48 * t69 - t58 * t91;
        const auto t93 = t70 * t90 + t72 * t92;
        const auto t94 = t50 * t51;
        const auto t95 = t65 * t94;
        const auto t96 = t70 * t89;
        const auto t97 = t72 * t91;
        const auto t98 = t75 * t94;
        const auto t99 = t16 * t91 - t20 * t48;
        const auto t100 = t30 * t89 - t41 * t47 + t99;
        const auto t101 = t73 * t94;
        const auto t102 = t73 * t78;
        const auto t103 = -t100 * t101 - t102 * t93 + t44 * t91 + t45 * t89
            + t93 * t95
            + t98 * (-t44 * t96 - t45 * t97 + 2 * t66 * t92 + 2 * t71 * t90);
        const auto t104 = t20 * t44 + t64 * t89;
        const auto t105 = t36 * t91 + t45 * t69;
        const auto t106 = t104 * t70 + t105 * t72;
        const auto t107 = t20 * t45 + t36 * t89;
        const auto t108 = t107 + t41 * t44 + t9 * t91;
        const auto t109 = -t101 * t108 + t102 * t106 - t106 * t95 + t47 * t91
            + t48 * t89
            - t98 * (t104 * t86 + t105 * t83 + t47 * t96 + t48 * t97);
        const auto t110 = t11 * t16 + t2 * t64;
        const auto t111 = t24 * t58 + t32 * t36;
        const auto t112 = t110 * t70 + t111 * t72;
        const auto t113 = t16 * t24 + t2 * t36;
        const auto t114 = t11 * t30 + t113 + t32 * t9;
        const auto t115 = t11 * t44;
        const auto t116 = t2 * t47;
        const auto t117 = t24 * t45;
        const auto t118 = t32 * t48;
        const auto t119 = t2 * t48 + t32 * t47;
        const auto t120 = t11 * t45 + t24 * t44;
        const auto t121 = -t101 * t114 + t102 * t112 - t112 * t95 + t119 + t120
            - t98
                * (t110 * t86 + t111 * t83 + t70 * (t115 + t116)
                   + t72 * (t117 + t118));
        const auto t122 = t16 * t5 - t2 * t20;
        const auto t123 = t27 * t58 - t32 * t69;
        const auto t124 = t122 * t70 + t123 * t72;
        const auto t125 = t16 * t27 - t2 * t69 - t20 * t32 + t5 * t58;
        const auto t126 = t16 + t30;
        const auto t127 = t102 * t124 - t124 * t50 * t51 * t65
            - t125 * t50 * t51 * t73 + t126 + t27 * t44 + t45 * t5
            - t50 * t51 * t75
                * (t122 * t86 + t123 * t83 + t70 * (t16 + t44 * t5)
                   + t72 * (t27 * t45 + t30));
        const auto t128 = t11 * t20 + t5 * t64;
        const auto t129 = t24 * t69 + t27 * t36;
        const auto t130 = t128 * t70 + t129 * t72;
        const auto t131 = t11 * t69 + t20 * t24 + t27 * t64 + t36 * t5;
        const auto t132 = t102 * t130 - t130 * t50 * t51 * t65
            - t131 * t50 * t51 * t73 + t27 * t47 + t36 + t48 * t5
            - t50 * t51 * t75
                * (t128 * t86 + t129 * t83 + t70 * (t47 * t5 + t9)
                   + t72 * (t27 * t48 + t36))
            + t9;
        const auto t133 = std::pow(t70 * t72, -1.0 / 2.0);
        const auto t134 = t13 * t44;
        const auto t135 = t47 * t6;
        const auto t136 = t134 + t135;
        const auto t137 = 1.0 / t72;
        const auto t138 = t13 * t16;
        const auto t139 = t138 + t6 * t64;
        const auto t140 = t36 * t6;
        const auto t141 = t13 * t58 + t140;
        const auto t142 = 1.0 / t70;
        const auto t143 = -t45;
        const auto t144 = -t58;
        const auto t145 = t52 * t67 - t55 * t62;
        const auto t146 = t143 * t144 + t145 * t48;
        const auto t147 = -t11;
        const auto t148 = t147 * t52 - t56 * t61;
        const auto t149 = -t64;
        const auto t150 = t56 * t59;
        const auto t151 = t147 * t62;
        const auto t152 = -t151;
        const auto t153 = t150 + t152;
        const auto t154 =
            std::pow(t148, 2) + std::pow(t149, 2) + std::pow(t153, 2);
        const auto t155 = -t47;
        const auto t156 = t148 * t44 + t149 * t155;
        const auto t157 = -t69;
        const auto t158 =
            std::pow(t144, 2) + std::pow(t145, 2) + std::pow(t157, 2);
        const auto t159 = t146 * t154 + t156 * t158;
        const auto t160 = std::pow(t72, -2);
        const auto t161 = t139 * t160;
        const auto t162 = t156 * t75;
        const auto t163 = 2 * t162;
        const auto t164 = t142 * t75;
        const auto t165 = t159 * t164;
        const auto t166 = t13 * t45;
        const auto t167 = t48 * t6;
        const auto t168 = t166 + t167;
        const auto t169 = t0 * t16 - t20 * t6;
        const auto t170 = -t169;
        const auto t171 = t137 * t65;
        const auto t172 = t0 * t44 + t16;
        const auto t173 = -t172;
        const auto t174 = t137 * t75;
        const auto t175 = t0 * t58 - t6 * t69;
        const auto t176 = t142 * t175;
        const auto t177 = t137 * t159;
        const auto t178 = 2 * t170;
        const auto t179 = t160 * t178;
        const auto t180 = t160 * t165;
        const auto t181 = t0 * t45 + t30;
        const auto t182 = t0 * t47 + t9;
        const auto t183 = t13 * t20;
        const auto t184 = t0 * t64 + t183;
        const auto t185 = t160 * t184;
        const auto t186 = t0 * t36;
        const auto t187 = t13 * t69 + t186;
        const auto t188 = t142 * t177;
        const auto t189 = t0 * t48 + t36;
        const auto t190 = -t171 * t184 + t180 * t184 - t187 * t188 + t189;
        const auto t191 = 2 * t141;
        const auto t192 = std::pow(t70, -2);
        const auto t193 = t192 * t75;
        const auto t194 = t146 * t193;
        const auto t195 = t174 * t192;
        const auto t196 = t159 * t195;
        const auto t197 = -t175;
        const auto t198 = 2 * t197;
        const auto t199 = t142 * t65;
        const auto t200 = t137 * t169;
        const auto t201 = t142 * t200;
        const auto t202 = -t159 * t201 + t172 - t196 * t197 + t197 * t199;
        const auto t203 = -t189;
        const auto t204 = 2 * t187;
        const auto t205 = 2 * t92;
        const auto t206 = t205 * t87;
        const auto t207 = 2 * t90;
        const auto t208 = t207 * t84;
        const auto t209 = 2 * t89 * t91;
        const auto t210 = t100 * t94;
        const auto t211 = 2 * t93;
        const auto t212 = std::pow(t89, 2);
        const auto t213 = std::pow(t91, 2);
        const auto t214 = t93 * t94;
        const auto t215 = t78 * t93;
        const auto t216 = t106 * t210 + t106 * t215 - t108 * t214 + t44 * t48
            + t45 * t47
            - t98
                * (t104 * t205 + t105 * t207 + t44 * t47 * t70
                   + t45 * t48 * t72);
        const auto t217 = -t12;
        const auto t218 = t14 + t217;
        const auto t219 = -t25;
        const auto t220 = t219 + t28;
        const auto t221 = -t11 * t91 + t112 * t210 + t112 * t215 - t114 * t214
            + t126 - t24 * t89
            + t98
                * (-t110 * t205 - t111 * t207 + t70 * (t11 * t89 + t218)
                   + t72 * (t220 + t24 * t91));
        const auto t222 = t5 * t89;
        const auto t223 = t27 * t91;
        const auto t224 = t27 * t89 + t5 * t91;
        const auto t225 = t119 - t124 * t210 - t124 * t215 + t125 * t214 + t224
            + t98
                * (2 * t122 * t92 + 2 * t123 * t90 - t70 * (t116 + t222)
                   - t72 * (t118 + t223));
        const auto t226 = t20 + t41;
        const auto t227 = t11 * t48 + t130 * t210 + t130 * t215
            - t131 * t50 * t51 * t93 + t226 + t24 * t47
            - t50 * t51 * t75
                * (t128 * t205 + t129 * t207 + t70 * (t11 * t47 + t20)
                   + t72 * (t24 * t48 + t41));
        const auto t228 = t13 * t89;
        const auto t229 = -t89;
        const auto t230 = t148 * t229 + t153 * t47;
        const auto t231 = t230 * t75;
        const auto t232 = 2 * t231;
        const auto t233 = t13 * t91;
        const auto t234 = -t47 * t69 + t58 * t89 + t99;
        const auto t235 = t137 * t234;
        const auto t236 = t144 * t91 - t157 * t48;
        const auto t237 = t154 * t236 + t158 * t230;
        const auto t238 = t164 * t237;
        const auto t239 =
            -t137 * t141 * t142 * t237 + t139 * t235 + t161 * t238 - t233 + t30;
        const auto t240 = t0 * t89;
        const auto t241 = t135 + t240;
        const auto t242 = t137 * t237;
        const auto t243 = t160 * t170;
        const auto t244 = t0 * t91;
        const auto t245 = t167 + t244;
        const auto t246 = t13 * t47 + t20;
        const auto t247 = t142 * t242;
        const auto t248 = t13 * t48 + t41;
        const auto t249 = t184 * t235 + t185 * t238 - t187 * t247 + t248;
        const auto t250 = t142 * t234;
        const auto t251 = -t220 - t233;
        const auto t252 = t193 * t236;
        const auto t253 = t195 * t237;
        const auto t254 = t16 - t228;
        const auto t255 = -t248;
        const auto t256 = 2 * t105;
        const auto t257 = t256 * t87;
        const auto t258 = 2 * t104;
        const auto t259 = t258 * t84;
        const auto t260 = -t3 + t7;
        const auto t261 = -t33 + t34;
        const auto t262 = t106 * t78;
        const auto t263 = -t106 * t114 * t50 * t51 - t108 * t112 * t50 * t51
            + t112 * t262 + t2 * t91 + t260 + t261 + t32 * t89
            - t50 * t51 * t75
                * (t110 * t256 + t111 * t258 + t70 * (t2 * t89 + t260)
                   + t72 * (t261 + t32 * t91));
        const auto t264 = -t17;
        const auto t265 = t18 + t264;
        const auto t266 = -t38;
        const auto t267 = t266 + t39;
        const auto t268 = t108 * t94;
        const auto t269 = t106 * t94;
        const auto t270 = t124 * t262 - t124 * t268 - t125 * t269 - t2 * t45
            + t226 - t32 * t44
            + t98
                * (-t122 * t256 - t123 * t258 + t70 * (t2 * t44 + t265)
                   + t72 * (t267 + t32 * t45));
        const auto t271 = t120 + t130 * t262 - t130 * t268 - t131 * t269 + t224
            - t98
                * (t128 * t256 + t129 * t258 + t70 * (t115 + t222)
                   + t72 * (t117 + t223));
        const auto t272 = t6 * t89;
        const auto t273 = t260 + t272;
        const auto t274 = t149 * t89 - t153 * t44;
        const auto t275 = t274 * t75;
        const auto t276 = 2 * t275;
        const auto t277 = t6 * t91;
        const auto t278 = t107 + t44 * t69 + t64 * t91;
        const auto t279 = t137 * t278;
        const auto t280 = -t145 * t91 + t157 * t45;
        const auto t281 = t154 * t280 + t158 * t274;
        const auto t282 = t164 * t281;
        const auto t283 =
            -t137 * t141 * t142 * t281 + t139 * t279 + t161 * t282 - t277 + t36;
        const auto t284 = t45 * t6;
        const auto t285 = t44 * t6;
        const auto t286 = -t265 - t285;
        const auto t287 = t137 * t281;
        const auto t288 = t134 + t240;
        const auto t289 = t142 * t287;
        const auto t290 = t142 * t278;
        const auto t291 = -t261 - t277;
        const auto t292 = t193 * t280;
        const auto t293 = t195 * t281;
        const auto t294 = t267 + t284;
        const auto t295 = t197 * t290 + t197 * t293 + t20 + t201 * t281 - t285;
        const auto t296 = t166 + t244;
        const auto t297 = 2 * t111;
        const auto t298 = t297 * t87;
        const auto t299 = 2 * t110;
        const auto t300 = t299 * t84;
        const auto t301 = 2 * t11 * t24;
        const auto t302 = 2 * t2 * t32;
        const auto t303 = std::pow(t11, 2);
        const auto t304 = std::pow(t2, 2);
        const auto t305 = std::pow(t24, 2);
        const auto t306 = std::pow(t32, 2);
        const auto t307 = t114 * t94;
        const auto t308 = t112 * t94;
        const auto t309 = t5 * t70;
        const auto t310 = t27 * t72;
        const auto t311 = t112 * t78;
        const auto t312 = t11 * t27 - t124 * t307 + t124 * t311 - t125 * t308
            + t24 * t5
            - t98 * (t11 * t309 + t122 * t297 + t123 * t299 + t24 * t310);
        const auto t313 = -t130 * t307 + t130 * t311 - t131 * t308 + t2 * t27
            + t32 * t5
            - t98 * (t128 * t297 + t129 * t299 + t2 * t309 + t310 * t32);
        const auto t314 = -t32 * t6;
        const auto t315 = t11 * t13;
        const auto t316 = t2 * t6;
        const auto t317 = t315 + t316;
        const auto t318 = t11 * t58 + t113 + t32 * t64;
        const auto t319 = t137 * t318;
        const auto t320 = t147 * t148 + t149 * t2;
        const auto t321 = t320 * t75;
        const auto t322 = 2 * t321;
        const auto t323 = t144 * t24 + t145 * t67;
        const auto t324 = t154 * t323 + t158 * t320;
        const auto t325 = t164 * t324;
        const auto t326 = 2 * t12;
        const auto t327 = t15 + t326;
        const auto t328 = t137 * t324;
        const auto t329 = 2 * t25;
        const auto t330 = t29 + t329;
        const auto t331 = 2 * t3;
        const auto t332 = t331 + t8;
        const auto t333 = t142 * t328;
        const auto t334 = 2 * t33;
        const auto t335 = -t334 + t34;
        const auto t336 = -t2 * t6;
        const auto t337 = t13 * t24;
        const auto t338 = t32 * t6;
        const auto t339 = t337 + t338;
        const auto t340 = t142 * t318;
        const auto t341 = t193 * t323;
        const auto t342 = t195 * t324;
        const auto t343 = t334 + t35;
        const auto t344 = -t331 + t7;
        const auto t345 = 2 * t123;
        const auto t346 = t345 * t87;
        const auto t347 = 2 * t73;
        const auto t348 = t122 * t84;
        const auto t349 = 2 * t106;
        const auto t350 = 2 * t112;
        const auto t351 = 2 * t27 * t5;
        const auto t352 = std::pow(t5, 2);
        const auto t353 = std::pow(t27, 2);
        const auto t354 = 2 * t348;
        const auto t355 = t124 * t94;
        const auto t356 = t130 * t94;
        const auto t357 = t11 * t32 - t124 * t130 * t78 + t125 * t356
            + t131 * t355 + t2 * t24
            + t98
                * (-t11 * t2 * t70 + 2 * t122 * t129 + 2 * t123 * t128
                   - t24 * t32 * t72);
        const auto t358 = 2 * t28;
        const auto t359 = 2 * t14;
        const auto t360 = t12 - t359;
        const auto t361 = -t360;
        const auto t362 = -t125;
        const auto t363 = t137 * t362;
        const auto t364 = t148 * t5 + t153 * t59;
        const auto t365 = t364 * t75;
        const auto t366 = 2 * t365;
        const auto t367 = t144 * t55 + t157 * t32;
        const auto t368 = t154 * t367 + t158 * t364;
        const auto t369 = t164 * t368;
        const auto t370 = t0 * t5;
        const auto t371 = t316 + t370;
        const auto t372 = t137 * t368;
        const auto t373 = 2 * t17;
        const auto t374 = t19 + t373;
        const auto t375 = 2 * t38;
        const auto t376 = t142 * t372;
        const auto t377 = t184 * t363 + t185 * t369 - t187 * t376 - t375 + t39;
        const auto t378 = t25 - t358;
        const auto t379 = -t378;
        const auto t380 = t142 * t362;
        const auto t381 = t193 * t367;
        const auto t382 = t195 * t368;
        const auto t383 = t0 * t27;
        const auto t384 = t338 + t383;
        const auto t385 = t375 + t40;
        const auto t386 = 2 * t129;
        const auto t387 = t386 * t87;
        const auto t388 = t128 * t84;
        const auto t389 = 2 * t388;
        const auto t390 = t131 * t137;
        const auto t391 = t3 - 2 * t7;
        const auto t392 = -t391;
        const auto t393 = t141 * t142;
        const auto t394 = t11 * t153 + t149 * t61;
        const auto t395 = t145 * t27 + t157 * t53;
        const auto t396 = t154 * t395 + t158 * t394;
        const auto t397 = t137 * t396;
        const auto t398 = t394 * t75;
        const auto t399 = 2 * t398;
        const auto t400 = t164 * t396;
        const auto t401 = t33 - 2 * t34;
        const auto t402 = 2 * t18;
        const auto t403 = -t17 + t402;
        const auto t404 = 2 * t39;
        const auto t405 = t266 + t404;
        const auto t406 = t315 + t370;
        const auto t407 = t142 * t397;
        const auto t408 = t337 + t383;
        const auto t409 = t131 * t142;
        const auto t410 = -t401;
        const auto t411 = t193 * t395;
        const auto t412 = t195 * t396;
        const auto t413 =
            -t131 * t142 * t197 + t197 * t412 + t201 * t396 + t264 + t402;
        const auto t414 = t138 + t6 * t9;
        const auto t415 = t414 * t84;
        const auto t416 = t139 * t84;
        const auto t417 = t13 * t148 + t149 * t62;
        const auto t418 = 2 * t417;
        const auto t419 = t160 * t418;
        const auto t420 = t414 * t50;
        const auto t421 = t141 * t51;
        const auto t422 = t421 * t50;
        const auto t423 = 2 * t416;
        const auto t424 = std::pow(t139, 2);
        const auto t425 = t139 * t174;
        const auto t426 = std::pow(t13, 2);
        const auto t427 = std::pow(t6, 2);
        const auto t428 = t426 + t427;
        const auto t429 = -2 * t139 * t141 - t428 * t75;
        const auto t430 = t133 * t137;
        const auto t431 = t170 * t174;
        const auto t432 = t0 * t13;
        const auto t433 = t432 * t75;
        const auto t434 = t139 * t175 - t141 * t170 + t170 * t425 + t433;
        const auto t435 = t0 * t6;
        const auto t436 = t139 * t187 + t141 * t184 + t435 * t75;
        const auto t437 = -t137 * t139 * t184 * t75 + t436;
        const auto t438 = std::pow(t141, 2);
        const auto t439 =
            t133 * (t137 * t424 + t142 * t438 - t393 * t425 - t428);
        const auto t440 = t175 * t51;
        const auto t441 = t420 * t75;
        const auto t442 =
            t43 * (-t169 * t420 - t175 * t421 + t432 + t440 * t441);
        const auto t443 = t0 * t9 + t183;
        const auto t444 = t187 * t51;
        const auto t445 =
            t43 * (-t187 * t421 - t420 * t443 + t435 + t441 * t444);
        const auto t446 = t169 * t50;
        const auto t447 = t51 * (t0 * t30 - t41 * t6);
        const auto t448 = t447 * t50;
        const auto t449 = t169 * t84;
        const auto t450 = t178 * t84;
        const auto t451 = t148 * t52 + t153 * t6;
        const auto t452 = 2 * t451;
        const auto t453 = std::pow(t0, 2);
        const auto t454 = t427 + t453;
        const auto t455 = -t454 * t75;
        const auto t456 = t174 * t178;
        const auto t457 = t175 * t184;
        const auto t458 = t170 * t187;
        const auto t459 = t13 * t6;
        const auto t460 = t459 * t75;
        const auto t461 = -t460;
        const auto t462 = t174 * t184;
        const auto t463 = t13 * t30 + t140;
        const auto t464 = t463 * t51;
        const auto t465 = t464 * t75;
        const auto t466 =
            t43 * (-t139 * t446 + t432 + t446 * t465 - t447 * t463);
        const auto t467 = t142 * t174;
        const auto t468 =
            -t133 * (t170 * t197 * t467 + t170 * t200 + t176 * t197 + t454);
        const auto t469 = t169 * t187;
        const auto t470 = t43 * (t187 * t447 + t443 * t446 + t459 - t469 * t98);
        const auto t471 = t0 * t149 + t153 * t56;
        const auto t472 = 2 * t471;
        const auto t473 = t160 * t472;
        const auto t474 = t13 * t41 + t186;
        const auto t475 = 3 * t106;
        const auto t476 = t184 * t84;
        const auto t477 = t184 * t50;
        const auto t478 = t474 * t51;
        const auto t479 = t478 * t50;
        const auto t480 = 3 * t476;
        const auto t481 = std::pow(t184, 2);
        const auto t482 = t426 + t453;
        const auto t483 = t184 * t204 + t482 * t75;
        const auto t484 =
            t43 * (-t139 * t477 + t435 - t463 * t478 + t465 * t477);
        const auto t485 = t43 * (t175 * t478 + t184 * t446 - t457 * t98 + t459);
        const auto t486 = std::pow(t187, 2);
        const auto t487 =
            t133 * (t137 * t481 + t142 * t486 - t184 * t187 * t467 - t482);
        const auto t488 = t463 * t87;
        const auto t489 = t191 * t87;
        const auto t490 = t144 * t56 + t145 * t6;
        const auto t491 = t164 * t490;
        const auto t492 = t133 * t142;
        const auto t493 =
            -t139 * t197 + t141 * t164 * t197 + t141 * t169 + t433;
        const auto t494 = -t141 * t142 * t187 * t75 + t436;
        const auto t495 = t0 * t144 + t157 * t62;
        const auto t496 = 2 * t495;
        const auto t497 = t446 * t51;
        const auto t498 = t198 * t87;
        const auto t499 = t175 * t87;
        const auto t500 = t164 * t495;
        const auto t501 = t164 * t197;
        const auto t502 = t184 * t197;
        const auto t503 = t187 * t87;
        const auto t504 = 3 * t503;
        const auto t505 = t164 * (t13 * t157 + t145 * t52);
        dA[0] = t43
            * (-t46 - t49 + 2 * t50 * t51 * t65 * t73
               + t50 * t51 * t75
                   * (4 * t66 * t71 + t70 * (t79 + t80) + t72 * (t81 + t82))
               - std::pow(t73, 2) * t78 - t73 * t85 - t73 * t88);
        dA[1] = t43 * (t103 - t85 * t93 - t88 * t93);
        dA[2] = t43 * (t106 * t85 + t106 * t88 + t109);
        dA[3] = t43 * (t112 * t85 + t112 * t88 + t121);
        dA[4] = t43 * (-t124 * t85 - t124 * t88 - t127);
        dA[5] = t43 * (-t130 * t85 - t130 * t88 - t132);
        dA[6] = t133
            * (t136 * t137 * t75 + t137 * t139 * t65 + t137 * t141 * t142 * t159
               - t161 * t163 - t161 * t165 - t168);
        dA[7] = t133
            * (-t162 * t179 + t170 * t171 - t170 * t180 + t173 * t174
               - t176 * t177 + t181);
        dA[8] = t133 * (t163 * t185 - t174 * t182 + t190);
        dA[9] = t133
            * (-t136 + t137 * t139 * t142 * t159 + t141 * t142 * t65
               - t141 * t196 + t142 * t168 * t75 - t191 * t194);
        dA[10] = t133 * (-t164 * t181 - t194 * t198 + t202);
        dA[11] = t133
            * (t164 * t203 + t182 - t184 * t188 + t187 * t196 - t187 * t199
               + t194 * t204);
        dA[12] = t43 * (t103 - t206 * t73 - t208 * t73);
        dA[13] = t43
            * (-t206 * t93 - t208 * t93 - t209 - t210 * t211 - t49
               + t50 * t51 * t75
                   * (t70 * (t212 + t80) + t72 * (t213 + t82) + 4 * t90 * t92)
               - t78 * std::pow(t93, 2));
        dA[14] = t43 * (t106 * t206 + t106 * t208 + t216);
        dA[15] = t43 * (t112 * t206 + t112 * t208 + t221);
        dA[16] = t43 * (-t124 * t206 - t124 * t208 + t225);
        dA[17] = t43 * (-t130 * t206 - t130 * t208 - t227);
        dA[18] = t133 * (t137 * t75 * (-t218 - t228) - t161 * t232 - t239);
        dA[19] = t133
            * (t137 * t241 * t75 - t170 * t235 - t176 * t242 - t179 * t231
               - t238 * t243 - t245);
        dA[20] = t133 * (-t174 * t246 + t185 * t232 + t249);
        dA[21] = t133
            * (t137 * t139 * t142 * t237 - t141 * t250 - t141 * t253
               + t142 * t251 * t75 - t191 * t252 - t254);
        dA[22] = t133
            * (t142 * t245 * t75 - t197 * t250 - t197 * t253 - t198 * t252
               - t201 * t237 - t241);
        dA[23] = t133
            * (t164 * t255 - t184 * t247 + t187 * t250 + t187 * t253
               + t204 * t252 + t246);
        dA[24] = t43 * (t109 + t257 * t73 + t259 * t73);
        dA[25] = t43 * (t216 + t257 * t93 + t259 * t93);
        dA[26] = t43
            * (-std::pow(t106, 2) * t78 + 2 * t106 * t108 * t50 * t51
               - t106 * t257 - t106 * t259 - t209 - t46
               + t50 * t51 * t75
                   * (4 * t104 * t105 + t70 * (t212 + t79)
                      + t72 * (t213 + t81)));
        dA[27] = t43 * (-t112 * t257 - t112 * t259 - t263);
        dA[28] = t43 * (t124 * t257 + t124 * t259 + t270);
        dA[29] = t43 * (t130 * t257 + t130 * t259 + t271);
        dA[30] = t133 * (-t137 * t273 * t75 - t161 * t276 - t283);
        dA[31] = t133
            * (t137 * t286 * t75 - t170 * t279 - t176 * t287 - t179 * t275
               - t243 * t282 + t284 - t41);
        dA[32] = t133
            * (-t166 + t174 * t288 + t184 * t279 + t185 * t276 + t185 * t282
               - t187 * t289 - t244);
        dA[33] = t133
            * (t137 * t139 * t142 * t281 - t141 * t290 - t141 * t293
               + t142 * t291 * t75 - t191 * t292 + t272 - t9);
        dA[34] = t133 * (-t142 * t294 * t75 - t198 * t292 - t295);
        dA[35] = t133
            * (-t134 + t164 * t296 - t184 * t289 + t187 * t290 + t187 * t293
               + t204 * t292 - t240);
        dA[36] = t43 * (t121 + t298 * t73 + t300 * t73);
        dA[37] = t43 * (t221 + t298 * t93 + t300 * t93);
        dA[38] = t43 * (-t106 * t298 - t106 * t300 - t263);
        dA[39] = t43
            * (-std::pow(t112, 2) * t78 + 2 * t112 * t114 * t50 * t51
               - t112 * t298 - t112 * t300 - t301 - t302
               + t50 * t51 * t75
                   * (4 * t110 * t111 + t70 * (t303 + t304)
                      + t72 * (t305 + t306)));
        dA[40] = t43 * (t124 * t298 + t124 * t300 + t312);
        dA[41] = t43 * (t130 * t298 + t130 * t300 + t313);
        dA[42] = t133
            * (t13 * t24 + t137 * t141 * t142 * t324 - t139 * t319 - t161 * t322
               - t161 * t325 - t174 * t317 - t314);
        dA[43] = t133
            * (t137 * t327 * t75 - t170 * t319 - t176 * t328 - t179 * t321
               - t243 * t325 - t330);
        dA[44] = t133
            * (t174 * t332 + t184 * t319 + t185 * t322 + t185 * t325
               - t187 * t333 + t335);
        dA[45] = t133
            * (t11 * t13 + t137 * t139 * t142 * t324 - t141 * t340 - t141 * t342
               - t164 * t339 - t191 * t341 - t336);
        dA[46] = t133
            * (t142 * t330 * t75 - t197 * t340 - t197 * t342 - t198 * t341
               - t201 * t324 - t327);
        dA[47] = t133
            * (t164 * t343 - t184 * t333 + t187 * t340 + t187 * t342
               + t204 * t341 + t344);
        dA[48] = t43 * (-t127 - t346 * t73 - t347 * t348);
        dA[49] = t43 * (-t211 * t348 + t225 - t346 * t93);
        dA[50] = t43 * (t106 * t346 + t270 + t348 * t349);
        dA[51] = t43 * (t112 * t346 + t312 + t348 * t350);
        dA[52] = t43
            * (-std::pow(t124, 2) * t78 + 2 * t124 * t125 * t50 * t51
               - t124 * t346 - t124 * t354 - t302 - t351
               + t50 * t51 * t75
                   * (4 * t122 * t123 + t70 * (t304 + t352)
                      + t72 * (t306 + t353)));
        dA[53] = t43 * (-t130 * t346 - t130 * t354 + t357);
        dA[54] = t133
            * (t137 * t141 * t142 * t368 + t137 * t361 * t75 - t139 * t363
               - t161 * t366 - t161 * t369 - t219 - t358);
        dA[55] = t133
            * (t0 * t27 - t170 * t363 - t174 * t371 - t176 * t372 - t179 * t365
               - t243 * t369 - t314);
        dA[56] = t133 * (t174 * t374 + t185 * t366 + t377);
        dA[57] = t133
            * (t137 * t139 * t142 * t368 - t141 * t380 - t141 * t382
               + t142 * t379 * t75 - t191 * t381 - t217 - t359);
        dA[58] = t133
            * (t0 * t5 - t164 * t384 - t197 * t380 - t197 * t382 - t198 * t381
               - t201 * t368 - t336);
        dA[59] = t133
            * (t164 * t385 + t18 - t184 * t376 + t187 * t380 + t187 * t382
               + t204 * t381 - t373);
        dA[60] = t43 * (-t132 - t347 * t388 - t387 * t73);
        dA[61] = t43 * (-t211 * t388 - t227 - t387 * t93);
        dA[62] = t43 * (t106 * t387 + t271 + t349 * t388);
        dA[63] = t43 * (t112 * t387 + t313 + t350 * t388);
        dA[64] = t43 * (-t124 * t387 - t124 * t389 + t357);
        dA[65] = t43
            * (-std::pow(t130, 2) * t78 + 2 * t130 * t131 * t50 * t51
               - t130 * t387 - t130 * t389 - t301 - t351
               + t50 * t51 * t75
                   * (4 * t128 * t129 + t70 * (t303 + t352)
                      + t72 * (t305 + t353)));
        dA[66] = t133
            * (t139 * t390 - t161 * t399 - t161 * t400 + t174 * t392
               + t393 * t397 + t401);
        dA[67] = t133
            * (t131 * t137 * t170 + t137 * t403 * t75 - t176 * t397
               - t179 * t398 - t243 * t400 - t405);
        dA[68] = t133
            * (-t174 * t406 - t184 * t390 + t185 * t399 + t185 * t400
               - t187 * t407 + t408);
        dA[69] = t133
            * (t139 * t407 + t141 * t409 - t141 * t412 + t164 * t410
               - t191 * t411 + t391);
        dA[70] = t133 * (t142 * t75 * (-t38 + t404) - t198 * t411 - t413);
        dA[71] = t133
            * (-t164 * t408 - t184 * t407 - t187 * t409 + t187 * t412
               + t204 * t411 + t406);
        dA[72] = t43
            * (t141 * t50 * t51 * t73 - t168 - t347 * t416 + t414 * t50 * t65
               - t415 * t73 + t50 * t51 * t75 * (t136 * t70 + t139 * t86));
        dA[73] = t133
            * (t137 * t142 * t75 * (t158 * (t13 * t229 + t148) + t236 * t418)
               - t238 * t419 - t239);
        dA[74] = t133
            * (t137 * t142 * t75 * (t158 * (t62 * t89 + t64) + t280 * t418)
               - t282 * t419 - t283);
        dA[75] = t43
            * (t112 * t415 - t112 * t422 - t114 * t420 + t339 + t350 * t416
               - t98 * (t139 * t297 + t317 * t70));
        dA[76] = t43
            * (-t124 * t415 + t124 * t422 - t124 * t423 + t125 * t420 + t378
               + t98 * (t139 * t345 + t361 * t70));
        dA[77] = t43
            * (-t130 * t415 + t130 * t422 - t130 * t423 + t131 * t420 + t401
               + t98 * (t139 * t386 + t392 * t70));
        dA[78] = t430 * (-t137 * t424 * t75 - t418 * t425 - t429);
        dA[79] = t430 * (-t418 * t431 - t434);
        dA[80] = t430 * (2 * t137 * t184 * t417 * t75 - t437);
        dA[81] = t439;
        dA[82] = t442;
        dA[83] = t445;
        dA[84] = t43
            * (t181 - t446 * t65 - t448 * t73 + t449 * t73 - t450 * t73
               + t98 * (t170 * t86 + t173 * t70));
        dA[85] = t43
            * (t100 * t169 * t50 + t169 * t51 * t75 * t76 * t93 - t245
               - t448 * t93 - t450 * t93
               + t50 * t51 * t75 * (t170 * t205 + t241 * t70));
        dA[86] = t43
            * (t106 * t448 - t106 * t449 + t106 * t450 + t108 * t446 + t294
               + t98 * (-t170 * t256 + t286 * t70));
        dA[87] = t43
            * (t112 * t448 - t112 * t449 + t112 * t450 + t114 * t446 + t28
               - t329 + t98 * (-t170 * t297 + t327 * t70));
        dA[88] = t43
            * (-t124 * t448 + t124 * t449 - t124 * t450 - t125 * t446 + t384
               + t98 * (2 * t123 * t170 - t371 * t70));
        dA[89] = t43
            * (t130 * t169 * t51 * t75 * t76 - t130 * t448 - t130 * t450
               - t131 * t446 - t405
               + t50 * t51 * t75 * (t170 * t386 + t403 * t70));
        dA[90] = t430 * (-t425 * t452 - t434);
        dA[91] = t430
            * (-std::pow(t170, 2) * t174 - t175 * t178 - t451 * t456 - t455);
        dA[92] = t430 * (t184 * t431 + t452 * t462 + t457 - t458 + t461);
        dA[93] = t466;
        dA[94] = t468;
        dA[95] = t470;
        dA[96] = t133
            * (-t180 * t472 + t190
               + t467 * (t146 * t472 + t158 * (t0 * t155 - t60 + t63)));
        dA[97] = t133
            * (-t238 * t473 + t249
               + t467 * (t158 * (-t150 + t151 + t47 * t56) + t236 * t472));
        dA[98] = t43
            * (t106 * t474 * t50 * t51 + t108 * t184 * t50 - t296 - t475 * t476
               + t50 * t51 * t75 * (t184 * t256 + t288 * t70));
        dA[99] = t43
            * (t112 * t479 - t112 * t480 + t114 * t477 + t335
               + t98 * (t184 * t297 + t332 * t70));
        dA[100] = t133
            * (-t369 * t473 + t377
               + t467 * (t158 * (2 * t150 + t152) + t367 * t472));
        dA[101] = t43
            * (-t130 * t479 + t130 * t480 - t131 * t477 + t408
               - t98 * (t184 * t386 + t406 * t70));
        dA[102] = t430 * (-t425 * t472 - t437);
        dA[103] = t430
            * (t137 * t170 * t184 * t75 + t175 * t184 - t456 * t471 - t458
               - t460);
        dA[104] = t430 * (-t137 * t481 * t75 + t462 * t472 + t483);
        dA[105] = t484;
        dA[106] = t485;
        dA[107] = t487;
        dA[108] = t43
            * (-t136 + t139 * t50 * t51 * t73 + t463 * t51 * t65 - t488 * t73
               - t489 * t73 + t50 * t51 * t75 * (t141 * t83 + t168 * t72));
        dA[109] = t43
            * (-t100 * t464 + t139 * t50 * t51 * t93 - t254 - t488 * t93
               - t489 * t93 + t50 * t51 * t75 * (t141 * t207 + t251 * t72));
        dA[110] = t43
            * (t106 * t488 + t106 * t489 - t108 * t464 - t139 * t269 + t273
               + t98 * (-t141 * t258 + t291 * t72));
        dA[111] = t43
            * (t112 * t488 + t112 * t489 - t114 * t464 - t139 * t308 + t317
               - t98 * (t141 * t299 + t339 * t72));
        dA[112] = t43
            * (-t124 * t488 - t124 * t489 + t125 * t464 + t139 * t355 + t360
               + t98 * (t122 * t191 + t379 * t72));
        dA[113] = t43
            * (-t130 * t488 - t130 * t489 + t131 * t464 + t139 * t356 + t391
               + t98 * (t128 * t191 + t410 * t72));
        dA[114] = t439;
        dA[115] = t466;
        dA[116] = t484;
        dA[117] = t492 * (-t142 * t438 * t75 - t191 * t491 - t429);
        dA[118] = t492 * (-t198 * t491 - t493);
        dA[119] = t492 * (2 * t142 * t187 * t490 * t75 - t494);
        dA[120] = t133
            * (-t196 * t496 + t202
               + t467 * (t154 * (t0 * t143 - t54 + t57) + t156 * t496));
        dA[121] = t43
            * (t100 * t175 * t51 + t175 * t50 * t74 * t75 * t93 - t241
               - t497 * t93 - t498 * t93
               + t50 * t51 * t75 * (t197 * t207 + t245 * t72));
        dA[122] = t133
            * (t137 * t142 * t75 * (t154 * (t45 * t62 + t69) + t274 * t496)
               - t293 * t496 - t295);
        dA[123] = t43
            * (t112 * t497 + t112 * t498 - t112 * t499 + t114 * t440 + t14
               - t326 + t98 * (-t197 * t299 + t330 * t72));
        dA[124] = t43
            * (-t124 * t497 - t124 * t498 + t124 * t499 - t125 * t440 + t371
               + t98 * (2 * t122 * t197 - t384 * t72));
        dA[125] = t133
            * (t137 * t142 * t75 * (t154 * (2 * t53 * t62 - t68) + t394 * t496)
               - t412 * t496 - t413);
        dA[126] = t442;
        dA[127] = t468;
        dA[128] = t485;
        dA[129] = t492 * (-t191 * t500 - t493);
        dA[130] = t492
            * (-t164 * std::pow(t197, 2) - t169 * t198 - t455 - t496 * t501);
        dA[131] = t492 * (t187 * t501 + t204 * t500 + t461 + t469 - t502);
        dA[132] = t43
            * (-t101 * t443 + t182 - t444 * t65 + t504 * t73
               + t98 * (-t187 * t83 + t203 * t72));
        dA[133] = t43
            * (t100 * t444 - t214 * t443 + t246 + t504 * t93
               + t98 * (-t187 * t207 + t255 * t72));
        dA[134] = t43
            * (t106 * t443 * t50 * t51 + t108 * t187 * t51 - t288 - t475 * t503
               + t50 * t51 * t75 * (t187 * t258 + t296 * t72));
        dA[135] = t43
            * (-t112 * t504 + t114 * t444 + t308 * t443 + t344
               + t98 * (t187 * t299 + t343 * t72));
        dA[136] = t43
            * (3 * t124 * t187 * t50 * t74 * t75 - t125 * t444 - t355 * t443
               - t374 + t50 * t51 * t75 * (-t122 * t204 + t385 * t72));
        dA[137] = t43
            * (t130 * t504 - t131 * t444 - t356 * t443 + t406
               - t98 * (t128 * t204 + t408 * t72));
        dA[138] = t445;
        dA[139] = t470;
        dA[140] = t487;
        dA[141] = t492 * (-t191 * t505 - t494);
        dA[142] = t492
            * (t142 * t187 * t197 * t75 + t169 * t187 - t198 * t505 - t460
               - t502);
        dA[143] = t492 * (-t142 * t486 * t75 + t204 * t505 + t483);
    }


    void normalized_mix_prod_gradient(
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
        double dA[12])
    {
        const auto t0 = t0_x - t1_x;
        const auto t1 = std::pow(t0, 2);
        const auto t2 = t0_y - t1_y;
        const auto t3 = std::pow(t2, 2);
        const auto t4 = t0_z - t1_z;
        const auto t5 = std::pow(t4, 2);
        const auto t6 = t1 + t3 + t5;
        const auto t7 = -t0;
        const auto t8 = -t2_y;
        const auto t9 = t0_y + t8;
        const auto t10 = -t9;
        const auto t11 = -t2_x;
        const auto t12 = t0_x + t11;
        const auto t13 = -t12;
        const auto t14 = -t2;
        const auto t15 = t10 * t7 - t13 * t14;
        const auto t16 = -t2_z;
        const auto t17 = t0_z + t16;
        const auto t18 = t0 * t17 - t12 * t4;
        const auto t19 = -t17;
        const auto t20 = -t4;
        const auto t21 = -t10 * t20 + t14 * t19;
        const auto t22 = std::pow(t15, 2) + std::pow(t18, 2) + std::pow(t21, 2);
        const auto t23 = -t3_z;
        const auto t24 = t0_z + t23;
        const auto t25 = -t24;
        const auto t26 = -t3_x;
        const auto t27 = t0_x + t26;
        const auto t28 = -t27;
        const auto t29 = -t20 * t28 + t25 * t7;
        const auto t30 = -t3_y;
        const auto t31 = t0_y + t30;
        const auto t32 = t0 * t31 - t2 * t27;
        const auto t33 = t2 * t24 - t31 * t4;
        const auto t34 = std::pow(t29, 2) + std::pow(t32, 2) + std::pow(t33, 2);
        const auto t35 = std::pow(t22 * t34 * t6, -1.0 / 2.0);
        const auto t36 = t18 * t32;
        const auto t37 = t15 * t29;
        const auto t38 = t1_y + t8;
        const auto t39 = t1_y + t30;
        const auto t40 = t16 + t1_z;
        const auto t41 = t1_z + t23;
        const auto t42 = -t38;
        const auto t43 = -t29;
        const auto t44 = -t41;
        const auto t45 = -t15;
        const auto t46 = -t13 * t20 + t19 * t7;
        const auto t47 = -t31;
        const auto t48 = -t14 * t28 + t47 * t7;
        const auto t49 = t14 * t25 - t20 * t47;
        const auto t50 = std::pow(t43, 2) + std::pow(t48, 2) + std::pow(t49, 2);
        const auto t51 = std::pow(t14, 2) + std::pow(t20, 2) + std::pow(t7, 2);
        const auto t52 = t50 * t51;
        const auto t53 = -t21;
        const auto t54 = std::pow(t45, 2) + std::pow(t46, 2) + std::pow(t53, 2);
        const auto t55 = t51 * t54;
        const auto t56 = t50 * t54;
        const auto t57 = 1.0 / t34;
        const auto t58 = t15 * t33;
        const auto t59 = t21 * t32;
        const auto t60 = t0 * (t43 * t45 - t46 * t48) - t2 * (-t58 + t59)
            + t4 * (-t43 * t53 + t46 * t49);
        const auto t61 = t60 / t22;
        const auto t62 = t57 * t61 / t6;
        const auto t63 = t11 + t1_x;
        const auto t64 = t1_x + t26;
        const auto t65 = -t64;
        const auto t66 = -t40;
        const auto t67 = -t58 + t59;
        const auto t68 = -t63;
        const auto t69 = -t39;
        const auto t70 = -t18 * t33 + t21 * t29;
        const auto t71 = t18 * t31;
        const auto t72 = t17 * t32;
        const auto t73 = t12 * t33 - t21 * t27;
        const auto t74 = t57 * t60;
        dA[0] = t35
            * (t0 * (-t39 * t46 - t40 * t48 + t42 * t43 + t44 * t45)
               - t2 * (t21 * t39 - t33 * t38) - t36 + t37
               - t4 * (t21 * t41 - t33 * t40)
               - t62
                   * (t0 * t56 + t52 * (t40 * t46 + t42 * t45)
                      + t55 * (t39 * t48 + t43 * t44)));
        dA[1] = t35
            * (-t0 * (-t18 * t64 + t29 * t63)
               + t2 * (-t41 * t45 + t48 * t66 - t49 * t63 + t53 * t65)
               - t4 * (-t18 * t41 + t29 * t40)
               - t62
                   * (t2 * t56 + t52 * (t45 * t63 + t53 * t66)
                      + t55 * (t41 * t49 + t48 * t65))
               - t67);
        dA[2] = t35
            * (-t0 * (t15 * t64 - t32 * t63) - t2 * (t15 * t39 - t32 * t38)
               + t4 * (-t38 * t43 + t46 * t69 + t49 * t68 - t53 * t64)
               - t62
                   * (t4 * t56 + t52 * (t38 * t53 + t46 * t68)
                      + t55 * (t43 * t64 + t49 * t69))
               - t70);
        dA[3] = t35
            * (-t0 * (t15 * t24 + t29 * t9 - t71 - t72)
               + t2 * (t47 * t53 - t49 * t9) + t36 - t37
               + t4 * (t19 * t49 - t24 * t53)
               - t62
                   * (t52 * (t19 * t46 + t45 * t9)
                      + t55 * (t24 * t43 + t47 * t48) + t56 * t7));
        dA[4] = t35
            * (t0 * (t13 * t43 - t27 * t46) - t2 * (t15 * t24 - t72 - t73)
               + t4 * (-t17 * t43 + t25 * t46)
               - t62
                   * (t14 * t56 + t52 * (t13 * t45 + t17 * t53)
                      + t55 * (t25 * t49 + t27 * t48))
               + t67);
        dA[5] = t35
            * (t0 * (-t12 * t48 + t28 * t45) + t2 * (t10 * t48 - t31 * t45)
               - t4 * (t29 * t9 - t71 - t73)
               - t62
                   * (t20 * t56 + t52 * (t10 * t53 + t12 * t46)
                      + t55 * (t28 * t43 + t31 * t49))
               + t70);
        dA[6] = t35
            * (t0 * (t14 * t43 - t4 * t48) + t3 * t33 + t33 * t5
               - t61 * (t15 * t2 + t18 * t4));
        dA[7] = -t35
            * (t1 * t29 + t2 * (t0 * t33 + t32 * t4) + t29 * t5
               + t61 * (-t0 * t15 + t21 * t4));
        dA[8] = t35
            * (t1 * t32 + t3 * t32 + t4 * (-t2 * t43 + t49 * t7)
               + t61 * (t0 * t18 + t2 * t21));
        dA[9] = -t35
            * (t0 * (-t15 * t4 + t18 * t2) + t21 * t3 + t21 * t5
               + t74 * (t2 * t32 + t29 * t4));
        dA[10] = t35
            * (t1 * t18 + t18 * t5 + t2 * (-t4 * t45 + t53 * t7)
               - t74 * (-t0 * t32 + t33 * t4));
        dA[11] = t35
            * (-t1 * t15 - t15 * t3 - t4 * (-t0 * t21 + t18 * t2)
               + t57 * t60 * (t0 * t29 + t2 * t33));
    }

    // dA is (144×1) flattened in column-major order
    void normalized_mix_prod_hessian(
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
        double dA[144])
    {
        const auto t0 = t0_x - t1_x;
        const auto t1 = std::pow(t0, 2);
        const auto t2 = t0_y - t1_y;
        const auto t3 = std::pow(t2, 2);
        const auto t4 = t0_z - t1_z;
        const auto t5 = std::pow(t4, 2);
        const auto t6 = t3 + t5;
        const auto t7 = t1 + t6;
        const auto t8 = -t2_y;
        const auto t9 = t0_y + t8;
        const auto t10 = t0 * t9;
        const auto t11 = -t2_x;
        const auto t12 = t0_x + t11;
        const auto t13 = t12 * t2;
        const auto t14 = -t13;
        const auto t15 = t10 + t14;
        const auto t16 = -t2_z;
        const auto t17 = t0_z + t16;
        const auto t18 = t0 * t17;
        const auto t19 = t12 * t4;
        const auto t20 = -t19;
        const auto t21 = t18 + t20;
        const auto t22 = std::pow(t21, 2);
        const auto t23 = t17 * t2;
        const auto t24 = t4 * t9;
        const auto t25 = -t24;
        const auto t26 = t23 + t25;
        const auto t27 = std::pow(t15, 2) + t22 + std::pow(t26, 2);
        const auto t28 = -t3_z;
        const auto t29 = t0_z + t28;
        const auto t30 = t0 * t29;
        const auto t31 = -t3_x;
        const auto t32 = t0_x + t31;
        const auto t33 = t32 * t4;
        const auto t34 = -t33;
        const auto t35 = t30 + t34;
        const auto t36 = -t3_y;
        const auto t37 = t0_y + t36;
        const auto t38 = t0 * t37;
        const auto t39 = t2 * t32;
        const auto t40 = -t39;
        const auto t41 = t38 + t40;
        const auto t42 = t2 * t29;
        const auto t43 = t37 * t4;
        const auto t44 = -t43;
        const auto t45 = t42 + t44;
        const auto t46 = std::pow(t41, 2) + std::pow(t45, 2);
        const auto t47 = std::pow(t35, 2) + t46;
        const auto t48 = std::pow(t27 * t47 * t7, -1.0 / 2.0);
        const auto t49 = t1_y + t8;
        const auto t50 = t1_z + t28;
        const auto t51 = t1_y + t36;
        const auto t52 = t16 + t1_z;
        const auto t53 = t49 * t50 - t51 * t52;
        const auto t54 = 2 * t0;
        const auto t55 = -t0;
        const auto t56 = -t9;
        const auto t57 = t55 * t56;
        const auto t58 = -t12;
        const auto t59 = -t2;
        const auto t60 = t58 * t59;
        const auto t61 = t57 - t60;
        const auto t62 = t21 * t52 + t49 * t61;
        const auto t63 = -t29;
        const auto t64 = t55 * t63;
        const auto t65 = -t32;
        const auto t66 = -t4;
        const auto t67 = t65 * t66;
        const auto t68 = t64 - t67;
        const auto t69 = t46 + std::pow(t68, 2);
        const auto t70 = t69 * t7;
        const auto t71 = t41 * t51 + t50 * t68;
        const auto t72 = -t17;
        const auto t73 = t59 * t72;
        const auto t74 = t56 * t66;
        const auto t75 = t73 - t74;
        const auto t76 = t22 + std::pow(t61, 2) + std::pow(t75, 2);
        const auto t77 = t7 * t76;
        const auto t78 = t69 * t76;
        const auto t79 = t0 * t78;
        const auto t80 = t62 * t70 + t71 * t77 + t79;
        const auto t81 = std::pow(t7, -2);
        const auto t82 = std::pow(t47, -2);
        const auto t83 = std::pow(t27, -2);
        const auto t84 = t61 * t68;
        const auto t85 = t21 * t41;
        const auto t86 = t45 * t61;
        const auto t87 = t41 * t75;
        const auto t88 = t86 - t87;
        const auto t89 = -t21 * t45;
        const auto t90 = t68 * t75 + t89;
        const auto t91 = t0 * (t84 - t85) + t2 * t88 - t4 * t90;
        const auto t92 = t83 * t91;
        const auto t93 = t81 * t82 * t92;
        const auto t94 = std::pow(t49, 2);
        const auto t95 = std::pow(t52, 2);
        const auto t96 = std::pow(t51, 2);
        const auto t97 = std::pow(t50, 2);
        const auto t98 = t62 * t7;
        const auto t99 = 4 * t0;
        const auto t100 = t62 * t69;
        const auto t101 = t71 * t76;
        const auto t102 = 1.0 / t47;
        const auto t103 = 1.0 / t27;
        const auto t104 = 1.0 / t7;
        const auto t105 = t103 * t104;
        const auto t106 = t102 * t105;
        const auto t107 = t106 * t91;
        const auto t108 = t45 * t49;
        const auto t109 = t51 * t75;
        const auto t110 = t108 - t109;
        const auto t111 = t45 * t52;
        const auto t112 = t50 * t75;
        const auto t113 = t111 - t112;
        const auto t114 = -t49;
        const auto t115 = -t68;
        const auto t116 = -t50;
        const auto t117 = -t61;
        const auto t118 = t55 * t72;
        const auto t119 = t58 * t66;
        const auto t120 = t118 - t119;
        const auto t121 = -t37;
        const auto t122 = t121 * t55;
        const auto t123 = t59 * t65;
        const auto t124 = t122 - t123;
        const auto t125 = -t84 + t85;
        const auto t126 =
            -t0 * (t114 * t115 + t116 * t117 - t120 * t51 - t124 * t52)
            - t110 * t2 - t113 * t4 + t125;
        const auto t127 = t106 * t80;
        const auto t128 = 2 * t71;
        const auto t129 = t105 * t82 * t91;
        const auto t130 = t128 * t129;
        const auto t131 = t102 * t80;
        const auto t132 = t104 * t92;
        const auto t133 = t132 * t62;
        const auto t134 = 2 * t133;
        const auto t135 = t54 * t81;
        const auto t136 = t103 * t91;
        const auto t137 = t135 * t136;
        const auto t138 = t41 * t52;
        const auto t139 = -2 * t138 + 2 * t15 * t50;
        const auto t140 = -2 * t21 * t51 + 2 * t35 * t49;
        const auto t141 = t2 * t78;
        const auto t142 = t11 + t1_x;
        const auto t143 = -t142 * t61 + t52 * t75;
        const auto t144 = t1_x + t31;
        const auto t145 = -t144 * t41 + t45 * t50;
        const auto t146 = t141 + t143 * t70 + t145 * t77;
        const auto t147 = t102 * t146;
        const auto t148 = t144 * t21;
        const auto t149 = t142 * t50 - t144 * t52;
        const auto t150 = t142 * t68;
        const auto t151 = -t148 + t150;
        const auto t152 = t21 * t50;
        const auto t153 = -t152;
        const auto t154 = t153 + t52 * t68;
        const auto t155 = -t144;
        const auto t156 = -t75;
        const auto t157 = -t52;
        const auto t158 = t59 * t63;
        const auto t159 = t121 * t66;
        const auto t160 = t158 - t159;
        const auto t161 = -t86 + t87;
        const auto t162 = t0 * t151 + t154 * t4 + t161
            - t2 * (-t117 * t50 + t124 * t157 - t142 * t160 + t155 * t156);
        const auto t163 = t106 * t126;
        const auto t164 = 2 * t2;
        const auto t165 = t100 * t164;
        const auto t166 = t101 * t164;
        const auto t167 = 2 * t98;
        const auto t168 = t143 * t7;
        const auto t169 = t142 * t70;
        const auto t170 = t144 * t77;
        const auto t171 = t54 * t69;
        const auto t172 = t54 * t76;
        const auto t173 = t143 * t171 + t145 * t172;
        const auto t174 = t80 * t93;
        const auto t175 = -t0 * t149
            - t107
                * (t128 * t168 + t145 * t167 + t165 + t166 - t169 * t49
                   - t170 * t51 + t173)
            + t110 + t127 * t162 + t146 * t163 + t146 * t174 + t148 - t150
            + t2 * t53;
        const auto t176 = t142 * t21 + t49 * t75;
        const auto t177 = t144 * t68 + t45 * t51;
        const auto t178 = -t4 * t78;
        const auto t179 = -t176 * t70 - t177 * t77 - t178;
        const auto t180 = 2 * t179;
        const auto t181 = t102 * t180;
        const auto t182 = t102 * t137;
        const auto t183 = -t51;
        const auto t184 = t142 * t51 - t144 * t49;
        const auto t185 = t142 * t41;
        const auto t186 = -t144 * t61 + t185;
        const auto t187 = t41 * t49;
        const auto t188 = t187 - t51 * t61;
        const auto t189 = -t142;
        const auto t190 = -t0 * t186 - t188 * t2
            - t4 * (-t115 * t49 + t120 * t183 - t144 * t156 + t160 * t189)
            + t90;
        const auto t191 = t128 * t7;
        const auto t192 = 2 * t4;
        const auto t193 = -t100 * t192 - t101 * t192;
        const auto t194 = t171 * t176 + t172 * t177;
        const auto t195 = t0 * t184
            - t107
                * (-t167 * t177 - t169 * t52 - t170 * t50 - t176 * t191 - t193
                   - t194)
            + t113 + t127 * t190 + t163 * t179 + t174 * t179 + t186
            + t4 * (-t116 * t49 + t183 * t52);
        const auto t196 = t17 * t21 + t61 * t9;
        const auto t197 = t29 * t68 + t37 * t41;
        const auto t198 = t196 * t70 + t197 * t77 + t79;
        const auto t199 = t182 * t198;
        const auto t200 = t45 * t9;
        const auto t201 = t37 * t75;
        const auto t202 = t200 - t201;
        const auto t203 = t17 * t45;
        const auto t204 = t203 - t29 * t75;
        const auto t205 = t21 * t37;
        const auto t206 = -t205;
        const auto t207 = t17 * t41;
        const auto t208 = -t207;
        const auto t209 = t206 + t208;
        const auto t210 = -t0 * (t15 * t29 + t209 + t35 * t9) - t15 * t35
            - t2 * t202 - t204 * t4 + t85;
        const auto t211 = t49 * t9;
        const auto t212 = t17 * t52;
        const auto t213 = t37 * t51;
        const auto t214 = t29 * t50;
        const auto t215 = t102 * t134;
        const auto t216 = t29 * t49;
        const auto t217 = t17 * t51;
        const auto t218 = t216 - t217;
        const auto t219 = t50 * t9;
        const auto t220 = t37 * t52;
        const auto t221 = t219 - t220;
        const auto t222 = t68 * t9;
        const auto t223 = t29 * t61;
        const auto t224 = t209 + t222 + t223;
        const auto t225 = -t21 * t51 + t49 * t68;
        const auto t226 = t50 * t61;
        const auto t227 = t226 - t41 * t52;
        const auto t228 = t0 * (t218 + t221) + t224 + t225 + t227;
        const auto t229 = -t116 * t17 + t52 * t63;
        const auto t230 = t32 * t52;
        const auto t231 = t12 * t50;
        const auto t232 = -t30 + t33;
        const auto t233 = t21 + t232;
        const auto t234 = t0 * (t230 - t231 + t233);
        const auto t235 = -t42 + t43;
        const auto t236 = t235 + t26;
        const auto t237 = t2 * (t218 + t236);
        const auto t238 = t12 * t61 - t17 * t75;
        const auto t239 = -t29 * t45 + t32 * t41;
        const auto t240 = -t141 + t238 * t70 + t239 * t77;
        const auto t241 = t32 * t51;
        const auto t242 = t12 * t68;
        const auto t243 = t21 * t32;
        const auto t244 = t242 - t243;
        const auto t245 = t17 * t68;
        const auto t246 = t21 * t29;
        const auto t247 = t245 - t246;
        const auto t248 = t12 * t45;
        const auto t249 = t32 * t75;
        const auto t250 = t248 - t249;
        const auto t251 = t207 - t223 + t250;
        const auto t252 =
            t0 * t244 - t15 * t45 + t2 * t251 + t247 * t4 + t26 * t41;
        const auto t253 = t32 * t49;
        const auto t254 = t12 * t51;
        const auto t255 = -t254;
        const auto t256 = -t38 + t39;
        const auto t257 = t15 + t256;
        const auto t258 = -t37 * t49 + t51 * t9;
        const auto t259 = t12 * t21 + t75 * t9;
        const auto t260 = t32 * t68 + t37 * t45;
        const auto t261 = t178 + t259 * t70 + t260 * t77;
        const auto t262 = t12 * t52;
        const auto t263 = t12 * t41;
        const auto t264 = t32 * t61;
        const auto t265 = t263 - t264;
        const auto t266 = t41 * t9;
        const auto t267 = t37 * t61;
        const auto t268 = t205 - t222 + t250;
        const auto t269 =
            -t0 * t265 + t2 * (-t266 + t267) + t26 * t35 + t268 * t4 + t89;
        const auto t270 = -t111 - t263 - t4 * (t221 + t236);
        const auto t271 = std::pow(t7 * t78, -1.0 / 2.0);
        const auto t272 = t4 * t51;
        const auto t273 = -t272;
        const auto t274 = t2 * t49;
        const auto t275 = t4 * t52;
        const auto t276 = t274 + t275;
        const auto t277 = t0 * (t115 * t117 - t120 * t124) + t2 * t88
            + t4 * (-t115 * t156 + t120 * t160);
        const auto t278 = 1.0 / t76;
        const auto t279 = t277 * t278;
        const auto t280 = t2 * t61;
        const auto t281 = t21 * t4;
        const auto t282 = t280 + t281;
        const auto t283 = -t126;
        const auto t284 = t278 * t283;
        const auto t285 = 2 * t282;
        const auto t286 = t114 * t117 + t120 * t52;
        const auto t287 = std::pow(t76, -2);
        const auto t288 = t277 * t287;
        const auto t289 = t286 * t288;
        const auto t290 = t3 * t45 + t45 * t5;
        const auto t291 = t0 * (t115 * t59 - t124 * t4) + t290;
        const auto t292 = 1.0 / t69;
        const auto t293 = t291 * t292;
        const auto t294 =
            std::pow(t115, 2) + std::pow(t124, 2) + std::pow(t160, 2);
        const auto t295 = std::pow(t55, 2);
        const auto t296 = std::pow(t59, 2);
        const auto t297 = std::pow(t66, 2);
        const auto t298 = t295 + t296 + t297;
        const auto t299 = t294 * t298;
        const auto t300 = t115 * t116 + t124 * t51;
        const auto t301 =
            std::pow(t117, 2) + std::pow(t120, 2) + std::pow(t156, 2);
        const auto t302 = t298 * t301;
        const auto t303 = t294 * t301;
        const auto t304 = t0 * t303 + t286 * t299 + t300 * t302;
        const auto t305 = t104 * t304;
        const auto t306 = t278 * t305;
        const auto t307 = t277 * t292;
        const auto t308 = t287 * t307;
        const auto t309 = t305 * t308;
        const auto t310 = t2 * t68;
        const auto t311 = t4 * t41;
        const auto t312 = -t311;
        const auto t313 = t310 + t312;
        const auto t314 = t0 * t68;
        const auto t315 = t0 * t61;
        const auto t316 = t4 * t75;
        const auto t317 = -t316;
        const auto t318 = t315 + t317;
        const auto t319 = -t318;
        const auto t320 = t0 * t49;
        const auto t321 = -t15 - t320;
        const auto t322 = t0 * t45;
        const auto t323 = t311 + t322;
        const auto t324 = t1 * t68 + t2 * t323 + t5 * t68;
        const auto t325 = t5 * t50;
        const auto t326 = -t158 + t159;
        const auto t327 = -t2 * (t326 + t51 * t66) + t325;
        const auto t328 = t2 * t50;
        const auto t329 = t3 * t51;
        const auto t330 = -t329;
        const auto t331 = t0 * t41;
        const auto t332 = 2 * t331;
        const auto t333 = t0 * t52;
        const auto t334 = t0 * t21;
        const auto t335 = t2 * t75;
        const auto t336 = t334 + t335;
        const auto t337 = 2 * t336;
        const auto t338 = t1 * t41 + t3 * t41;
        const auto t339 = t338 + t4 * (-t115 * t2 + t160 * t55);
        const auto t340 = t292 * t339;
        const auto t341 = -t278 * t283 * t336 + t306 * t340 + t309 * t336;
        const auto t342 = t2 * t52;
        const auto t343 = t4 * t49;
        const auto t344 = t0 * (t342 - t343);
        const auto t345 = t2 * t51;
        const auto t346 = t4 * t50;
        const auto t347 = t345 + t346;
        const auto t348 = t2 * t41;
        const auto t349 = t4 * t68;
        const auto t350 = t348 + t349;
        const auto t351 = t283 * t292;
        const auto t352 = std::pow(t69, -2);
        const auto t353 = t2 * t21;
        const auto t354 = t4 * t61;
        const auto t355 = t353 - t354;
        const auto t356 = t0 * t355 + t3 * t75 + t5 * t75;
        const auto t357 = t5 * t52;
        const auto t358 = 2 * t334;
        const auto t359 = t0 * t51;
        const auto t360 = t331 - t4 * t45;
        const auto t361 = -t360;
        const auto t362 = t351 * t361;
        const auto t363 = 2 * t361;
        const auto t364 = t277 * t352;
        const auto t365 = t300 * t364;
        const auto t366 = t1 * t21 + t21 * t5;
        const auto t367 = t2 * (-t117 * t4 + t156 * t55) + t366;
        const auto t368 = t278 * t367;
        const auto t369 = t292 * t305;
        const auto t370 = t368 * t369;
        const auto t371 = t279 * t361;
        const auto t372 = t305 * t352;
        const auto t373 = t0 * t50;
        const auto t374 = -t35 - t373;
        const auto t375 = t2 * t45;
        const auto t376 = t314 + t375;
        const auto t377 = 2 * t376;
        const auto t378 = t0 * t75;
        const auto t379 = -t353;
        const auto t380 = t378 + t379;
        const auto t381 = t1 * t61 + t3 * t61 - t380 * t4;
        const auto t382 = t279 * t372;
        const auto t383 = t3 * t49;
        const auto t384 = t383 - t4 * (t52 * t59 + t75);
        const auto t385 = 2 * t145;
        const auto t386 = t129 * t385;
        const auto t387 = t132 * t143;
        const auto t388 = 2 * t387;
        const auto t389 = t136 * t81;
        const auto t390 = t164 * t389;
        const auto t391 = -t149;
        const auto t392 = std::pow(t142, 2);
        const auto t393 = std::pow(t144, 2);
        const auto t394 = 4 * t2;
        const auto t395 = t143 * t69;
        const auto t396 = t145 * t76;
        const auto t397 = t106 * t146;
        const auto t398 = t142 * t45;
        const auto t399 = 2 * t144 * t26 - 2 * t398;
        const auto t400 = t102 * t390;
        const auto t401 = t35 * t52;
        const auto t402 = t106 * t162;
        const auto t403 = t176 * t7;
        const auto t404 = 2 * t168;
        const auto t405 = -t192 * t395 - t192 * t396;
        const auto t406 = t164 * t69;
        const auto t407 = t164 * t76;
        const auto t408 = t176 * t406 + t177 * t407;
        const auto t409 = t146 * t93;
        const auto t410 = -t107
                * (-t177 * t404 - t385 * t403 - t405 - t408 - t49 * t52 * t70
                   - t50 * t51 * t77)
            - t15 * t51 + t152 + t179 * t402 + t179 * t409 + t184 * t2 + t187
            + t190 * t397 + t391 * t4 - t401;
        const auto t411 = t142 * t29;
        const auto t412 = t144 * t17;
        const auto t413 = t411 - t412;
        const auto t414 = -t219;
        const auto t415 = t102 * t388;
        const auto t416 = t142 * t9;
        const auto t417 = -t10 + t13;
        const auto t418 = t385 * t7;
        const auto t419 = t240 * t400;
        const auto t420 = -t230 + t231;
        const auto t421 = t2 * (-t413 - t420);
        const auto t422 = t144 * t75;
        const auto t423 = t12 * t142;
        const auto t424 = t144 * t32;
        const auto t425 = -t266;
        const auto t426 = t144 * t9;
        const auto t427 = t142 * t37;
        const auto t428 = -t427;
        const auto t429 = t37 * t50;
        const auto t430 = 2 * t375;
        const auto t431 = -t162;
        const auto t432 = t278 * t431;
        const auto t433 = t142 * t2;
        const auto t434 = -t417 - t433;
        const auto t435 = t117 * t142 + t156 * t157;
        const auto t436 = t288 * t435;
        const auto t437 = t124 * t155 + t160 * t50;
        const auto t438 = t2 * t303 + t299 * t435 + t302 * t437;
        const auto t439 = t104 * t438;
        const auto t440 = t278 * t439;
        const auto t441 = t308 * t439;
        const auto t442 = t144 * t4;
        const auto t443 = -t0 * (-t35 - t442) + t325;
        const auto t444 = t0 * t142;
        const auto t445 = t275 + t444;
        const auto t446 = t2 * (t373 - t442) + t323;
        const auto t447 = t1 * t144;
        const auto t448 = -t26 - t342;
        const auto t449 = -t118 + t119;
        const auto t450 = t144 * t2;
        const auto t451 = -t256 - t450;
        const auto t452 = t292 * t431;
        const auto t453 = t142 * t4;
        const auto t454 = -t453;
        const auto t455 = t0 * t144;
        const auto t456 = t346 + t455;
        const auto t457 = t364 * t437;
        const auto t458 = t292 * t439;
        const auto t459 = t352 * t439;
        const auto t460 = t354 + t378;
        const auto t461 = t1 * t142;
        const auto t462 = 2 * t280;
        const auto t463 = -t18 + t19;
        const auto t464 = -t328 - t45;
        const auto t465 = t279 * t459;
        const auto t466 = t376 * t465;
        const auto t467 = 2 * t177;
        const auto t468 = t129 * t467;
        const auto t469 = t132 * t176;
        const auto t470 = 2 * t469;
        const auto t471 = t192 * t389;
        const auto t472 = 4 * t4;
        const auto t473 = t106 * t190;
        const auto t474 = t102 * t471;
        const auto t475 = -t426 + t427;
        const auto t476 = t0 * (t257 + t475);
        const auto t477 = -t216;
        const auto t478 = t106 * t179;
        const auto t479 = t144 * t29;
        const auto t480 = t192 * t69;
        const auto t481 = t192 * t76;
        const auto t482 = t196 * t7;
        const auto t483 = 2 * t403;
        const auto t484 = t17 * t49;
        const auto t485 = -t23 + t24;
        const auto t486 = t467 * t7;
        const auto t487 = t102 * t470;
        const auto t488 = -t253 + t254;
        const auto t489 = -t0 * (t144 * t58 - t189 * t32) + t188
            + t2 * (t257 + t488) - t245 + t246 + t4 * (t233 - t411 + t412);
        const auto t490 = t261 * t474;
        const auto t491 = -t142 * t45 - t248 + t249 + t422;
        const auto t492 = t206 + t222 + t225 + t4 * (t475 + t488) + t491;
        const auto t493 = -t122 + t123;
        const auto t494 = -t453 - t463;
        const auto t495 = -t190;
        const auto t496 = t278 * t495;
        const auto t497 = t120 * t189 + t156 * t49;
        const auto t498 = t115 * t144 + t160 * t183;
        const auto t499 = t299 * t497 + t302 * t498 + t303 * t4;
        const auto t500 = t104 * t499;
        const auto t501 = t278 * t500;
        const auto t502 = -2 * t349;
        const auto t503 = -t343 - t485;
        const auto t504 = 2 * t319;
        const auto t505 = t288 * t497;
        const auto t506 = t292 * t324;
        const auto t507 = t308 * t500;
        const auto t508 = -t319 * t496 + t319 * t507 + t501 * t506;
        const auto t509 = -t450;
        const auto t510 = t274 + t444;
        const auto t511 = -t310;
        const auto t512 = t322 + t511;
        const auto t513 = -2 * t316;
        const auto t514 = t292 * t495;
        const auto t515 = -t232 - t442;
        const auto t516 = 2 * t350;
        const auto t517 = t364 * t498;
        const auto t518 = t278 * t356;
        const auto t519 = t292 * t500;
        const auto t520 = t279 * t352 * t500;
        const auto t521 = -t0 * (-t15 - t433) + t383;
        const auto t522 = -t235 - t272;
        const auto t523 = t345 + t455;
        const auto t524 = t353 - t378 + t4 * (t320 - t433);
        const auto t525 = t115 * t29 + t121 * t124;
        const auto t526 = 2 * t298;
        const auto t527 = t525 * t526;
        const auto t528 = t117 * t9 + t120 * t72;
        const auto t529 = t526 * t528;
        const auto t530 = 2 * t55;
        const auto t531 = t294 * t530;
        const auto t532 = t301 * t530;
        const auto t533 = t301 * t525;
        const auto t534 = t294 * t528;
        const auto t535 = -t303;
        const auto t536 = t104 * t292;
        const auto t537 = t279 * t536;
        const auto t538 = -t0 * t224 + t125 + t2 * (t121 * t156 - t160 * t9)
            + t4 * (-t156 * t29 + t160 * t72);
        const auto t539 = t278 * t538;
        const auto t540 = t299 * t528 + t302 * t525 + t303 * t55;
        const auto t541 = t536 * t540;
        const auto t542 = t304 * t81;
        const auto t543 = t279 * t292;
        const auto t544 = t542 * t543;
        const auto t545 = -t64 + t67;
        const auto t546 = t120 + t545;
        const auto t547 = t326 + t75;
        const auto t548 = -t57 + t60;
        const auto t549 = t352 * t540;
        const auto t550 = t288 * t549 * t81;
        const auto t551 = t135 * t543;
        const auto t552 = 2 * t528;
        const auto t553 = 2 * t525;
        const auto t554 = t17 * t37 - t29 * t9;
        const auto t555 = 2 * t205;
        const auto t556 = 2 * t207;
        const auto t557 = 2 * t106;
        const auto t558 = 2 * t198;
        const auto t559 = t102 * t132;
        const auto t560 = std::pow(t9, 2);
        const auto t561 = std::pow(t17, 2);
        const auto t562 = std::pow(t37, 2);
        const auto t563 = std::pow(t29, 2);
        const auto t564 = t117 * t58 + t156 * t17;
        const auto t565 = t124 * t32 + t160 * t63;
        const auto t566 = t299 * t564 + t302 * t565 + t303 * t59;
        const auto t567 = t104 * t308;
        const auto t568 = t566 * t567;
        const auto t569 = t104 * t279 * t352 * t566;
        const auto t570 = 2 * t59;
        const auto t571 = t0 * (t115 * t58 - t120 * t32) + t161 + t2 * t251
            + t4 * (-t115 * t17 + t120 * t63);
        const auto t572 = t278 * t571;
        const auto t573 = t536 * t539;
        const auto t574 = -t242 + t243;
        const auto t575 = t0 * (t29 * t58 - t32 * t72)
            + t2 * (t121 * t17 - t63 * t9) + t202
            + t4 * (-t17 * t29 + t63 * t72)
            - t537
                * (t121 * t302 * t32 + t299 * t58 * t9 + t527 * t564
                   + t529 * t565 + t531 * t564 + t532 * t565 + t533 * t570
                   + t534 * t570)
            - t541 * t572 + t550 * t566 - t566 * t573 + t574;
        const auto t576 = t115 * t65 + t160 * t37;
        const auto t577 = t12 * t120 + t156 * t56;
        const auto t578 = t299 * t577 + t302 * t576 + t303 * t66;
        const auto t579 = t567 * t578;
        const auto t580 = t352 * t578;
        const auto t581 = t104 * t279 * t580;
        const auto t582 = 2 * t66;
        const auto t583 = t0 * (t117 * t65 - t12 * t124)
            + t2 * (-t117 * t37 + t124 * t56) + t268 * t4 + t90;
        const auto t584 = t278 * t583;
        const auto t585 = t0 * (-t12 * t121 + t65 * t9)
            + t2 * (t121 * t56 - t37 * t9) + t204 + t265 - t4 * t554
            - t537
                * (t12 * t299 * t72 + t29 * t302 * t65 + t527 * t577
                   + t529 * t576 + t531 * t577 + t532 * t576 + t533 * t582
                   + t534 * t582)
            - t541 * t584 + t550 * t578 - t573 * t578;
        const auto t586 = t2 * t9;
        const auto t587 = t17 * t4;
        const auto t588 = t586 + t587;
        const auto t589 = t288 * t528;
        const auto t590 = t104 * t278;
        const auto t591 = t540 * t590;
        const auto t592 = t540 * t567;
        const auto t593 = t0 * t115;
        const auto t594 = 2 * t10 + t14;
        const auto t595 = 2 * t18 + t20;
        const auto t596 = t3 * t37 + t37 * t5;
        const auto t597 = -t0 * t26;
        const auto t598 = t2 * t37;
        const auto t599 = t29 * t4;
        const auto t600 = t598 + t599;
        const auto t601 = t292 * t538;
        const auto t602 = t364 * t525;
        const auto t603 = t104 * t279 * t549;
        const auto t604 = 2 * t38 + t40;
        const auto t605 = t17 * t3 + t17 * t5;
        const auto t606 = t0 * t117;
        const auto t607 = 2 * t30 + t34;
        const auto t608 = t526 * t564;
        const auto t609 = t526 * t565;
        const auto t610 = t294 * t564;
        const auto t611 = t301 * t565;
        const auto t612 = t294 * t570;
        const auto t613 = t301 * t570;
        const auto t614 = t536 * t566;
        const auto t615 = t543 * t81;
        const auto t616 = t164 * t615;
        const auto t617 = -t73 + t74;
        const auto t618 = t12 * t29 - t17 * t32;
        const auto t619 = std::pow(t12, 2);
        const auto t620 = std::pow(t32, 2);
        const auto t621 = 4 * t7;
        const auto t622 = 2 * t248 - 2 * t32 * t75;
        const auto t623 = -t12 * t37 + t32 * t9;
        const auto t624 = t536 * t578;
        const auto t625 = t267 + t425;
        const auto t626 = -t0 * (-t12 * t32 + t58 * t65) + t2 * t623 + t247
            - t277 * t287 * t352 * t566 * t578 * t81 + t4 * t618
            + t537
                * (t17 * t299 * t56 + t302 * t37 * t63 + t576 * t608
                   + t576 * t613 + t577 * t609 + t577 * t612 + t582 * t610
                   + t582 * t611)
            + t572 * t624 + t584 * t614 + t625;
        const auto t627 = t29 * t3;
        const auto t628 = -t10 + 2 * t12 * t2;
        const auto t629 = t566 * t590;
        const auto t630 = t29 * t5;
        const auto t631 = t1 * t29 + t630;
        const auto t632 = t0 * t12;
        const auto t633 = t587 + t632;
        const auto t634 = t288 * t564;
        const auto t635 = t124 * t2;
        const auto t636 = 2 * t23 + t25;
        const auto t637 = t295 * t72;
        const auto t638 = t156 * t2;
        const auto t639 = 2 * t2 * t32 - t38;
        const auto t640 = t292 * t571;
        const auto t641 = t350 * t640;
        const auto t642 = t0 * t32;
        const auto t643 = t599 + t642;
        const auto t644 = -t2 * (t17 * t55 - t4 * t58) + t460;
        const auto t645 = t12 * t3;
        const auto t646 = 2 * t42 + t44;
        const auto t647 = t1 * t12;
        const auto t648 = t12 * t5 + t647;
        const auto t649 = t493 + t61;
        const auto t650 = t526 * t576;
        const auto t651 = t526 * t577;
        const auto t652 = t301 * t576;
        const auto t653 = t294 * t577;
        const auto t654 = t294 * t582;
        const auto t655 = t301 * t582;
        const auto t656 = t288 * t580;
        const auto t657 = 2 * t577;
        const auto t658 = 2 * t576;
        const auto t659 = t192 * t615;
        const auto t660 = 2 * t12 * t4 - t18;
        const auto t661 = t160 * t4;
        const auto t662 = t578 * t590;
        const auto t663 = -t104 * t277 * t282 * t287 * t292 * t578 + t121 * t295
            + t282 * t584 + t293 * t662 - t296 * t37 + t661;
        const auto t664 = t3 * t32;
        const auto t665 = -t23 + 2 * t4 * t9;
        const auto t666 = t1 * t32;
        const auto t667 = t32 * t5 + t666;
        const auto t668 = t586 + t632;
        const auto t669 = t1 * t9;
        const auto t670 = t292 * t583;
        const auto t671 = -t30 + 2 * t32 * t4;
        const auto t672 = t3 * t9;
        const auto t673 = t5 * t9 + t672;
        const auto t674 = t120 * t4;
        const auto t675 = 2 * t37 * t4 - t42;
        const auto t676 = t598 + t642;
        const auto t677 = -t15 * t4 + t380;
        const auto t678 = t2 * t35;
        const auto t679 = t15 * t2 + t281;
        const auto t680 = t103 * t679;
        const auto t681 = t102 * (t0 * t313 + t290);
        const auto t682 = t105 * t681;
        const auto t683 = t285 * t7;
        const auto t684 = t171 * t282;
        const auto t685 = t132 * t679;
        const auto t686 = t132 * t285;
        const auto t687 = t282 * t406;
        const auto t688 = t102 * t685;
        const auto t689 = t285 * t559;
        const auto t690 = -t678;
        const auto t691 = t117 * t59;
        const auto t692 = t674 + t691;
        const auto t693 = -t277 * t6;
        const auto t694 = t279 * t692;
        const auto t695 = t271 * t278;
        const auto t696 = t0 * t277;
        const auto t697 = t2 * t696;
        const auto t698 = t279 * t282;
        const auto t699 = t282 * t324 - t291 * t319 + t319 * t698 + t697;
        const auto t700 = t282 * t339;
        const auto t701 = -t0 * t277 * t4;
        const auto t702 = t336 * t698;
        const auto t703 = t348 + t35 * t4;
        const auto t704 = t102 * t703;
        const auto t705 = t680 * t91;
        const auto t706 = t48 * (t356 * t680 - t681 * t703 + t704 * t705);
        const auto t707 = t2 * t460 + t366;
        const auto t708 = t102 * t360;
        const auto t709 =
            t48 * (t360 * t681 + t4 * t7 - t680 * t707 - t705 * t708);
        const auto t710 = t0 * t55;
        const auto t711 = t4 * t66;
        const auto t712 = t271
            * (t278 * t282 * t381 + t291 * t292 * t376 - t292 * t376 * t698
               - t59 * (-t296 + t710 + t711));
        const auto t713 = t103 * t318;
        const auto t714 = t504 * t7;
        const auto t715 = t171 * t319;
        const auto t716 = t132 * t318;
        const auto t717 = t319 * t406;
        const auto t718 = t102 * t716;
        const auto t719 = t115 * t66;
        const auto t720 = t156 * t66 + t606;
        const auto t721 = t102 * t324;
        const auto t722 = t105 * t721;
        const auto t723 = t504 * t559;
        const auto t724 = t279 * t720;
        const auto t725 = -t277 * (t1 + t5);
        const auto t726 = t279 * t319;
        const auto t727 =
            -t2 * t277 * t4 + t319 * t339 + t324 * t336 + t336 * t726;
        const auto t728 = t2 * t59;
        const auto t729 = t271
            * (t292 * t350 * t726 + t319 * t518 + t350 * t506
               - t66 * (-t297 + t710 + t728));
        const auto t730 = t713 * t91;
        const auto t731 = t48 * (-t360 * t721 + t707 * t713 + t708 * t730);
        const auto t732 = t102 * t376;
        const auto t733 =
            t48 * (t0 * t7 - t376 * t721 - t381 * t713 + t730 * t732);
        const auto t734 = t124 * t55;
        const auto t735 = t493 + t51 * t55;
        const auto t736 = t449 + t52 * t55;
        const auto t737 = t120 * t55;
        const auto t738 = t638 + t737;
        const auto t739 = -t348;
        const auto t740 = t103 * t336;
        const auto t741 = t338 - t4 * t512;
        const auto t742 = t102 * t741;
        const auto t743 = t105 * t742;
        const auto t744 = t336 * t406;
        const auto t745 = t337 * t7;
        const auto t746 = 3 * t336;
        const auto t747 = -t336 * t480;
        const auto t748 = t559 * t746;
        const auto t749 = t4 * t696;
        const auto t750 = t279 * t738;
        const auto t751 = t277 * (t1 + t3);
        const auto t752 = t740 * t91;
        const auto t753 =
            t48 * (t2 * t7 - t356 * t740 - t703 * t742 - t704 * t752);
        const auto t754 = t271
            * (t278 * t336 * t367 - t292 * t336 * t371 - t340 * t361
               - t55 * (-t295 + t711 + t728));
        const auto t755 = t48 * (t376 * t742 - t381 * t740 + t732 * t752);
        const auto t756 = t172 * t350;
        const auto t757 = t129 * t703;
        const auto t758 = t129 * t516;
        const auto t759 = t354 + t379;
        const auto t760 = t350 * t481;
        const auto t761 = t106 * t356;
        const auto t762 = t635 + t719;
        const auto t763 = t259 * t7;
        const auto t764 = t307 * t762;
        const auto t765 = t271 * t292;
        const auto t766 = t307 * t350;
        const auto t767 = -t350 * t367 + t356 * t361 + t361 * t766 + t697;
        const auto t768 = t356 * t376;
        const auto t769 = t376 * t766;
        const auto t770 = t661 + t734;
        const auto t771 = t361 * t407;
        const auto t772 = t129 * t363;
        const auto t773 = t129 * t360;
        const auto t774 = -t281;
        const auto t775 = t361 * t481;
        const auto t776 = t106 * t707;
        const auto t777 = t307 * t770;
        const auto t778 =
            t2 * t277 * t4 - t307 * t361 * t376 + t361 * t381 + t367 * t376;
        const auto t779 = t172 * t376;
        const auto t780 = 3 * t129 * t376;
        const auto t781 = t160 * t59 + t593;
        const auto t782 = -t376 * t481;
        const auto t783 = t307 * t781;
        dA[0] = t48
            * (-t107
                   * (t100 * t99 + t101 * t99 + t70 * (t94 + t95)
                      + 4 * t71 * t98 + t77 * (t96 + t97) + t78)
               + 2 * t126 * t127 + t130 * t80 + t131 * t134 + t131 * t137 + t139
               + t140 + t53 * t54 + std::pow(t80, 2) * t93);
        dA[1] = t48 * (t130 * t146 + t134 * t147 + t137 * t147 + t175);
        dA[2] = t48 * (t130 * t179 + t133 * t181 + t179 * t182 + t195);
        dA[3] = t48
            * (t102 * t103 * t104 * t91
                   * (t100 * t54 + t101 * t54 + t167 * t197 + t171 * t196
                      + t172 * t197 + t191 * t196 + t70 * (t211 + t212)
                      + t77 * (t213 + t214) + t78)
               - t127 * t210 - t130 * t198 - t163 * t198 - t174 * t198
               - t198 * t215 - t199 - t228);
        dA[4] = t48
            * (-t107
                   * (-t165 - t166 + t167 * t239 + t171 * t238 + t172 * t239
                      + t191 * t238 + t70 * (t12 * t49 + t15)
                      + t77 * (t241 + t41))
               - t108 + t109 - t127 * t252 + t130 * t240 + t163 * t240
               + t174 * t240 + t182 * t240 + t215 * t240 + t229 * t4 - t234
               - t237 + t244);
        dA[5] = t48
            * (t0 * (t253 + t255 + t257)
               - t107
                   * (t167 * t260 + t171 * t259 + t172 * t260 + t191 * t259
                      + t193 + t70 * (t21 + t262) + t77 * (t32 * t50 + t35))
               - t127 * t269 + t130 * t261 + t15 * t32 + t163 * t261
               + t174 * t261 + t182 * t261 - t2 * t258 + t215 * t261 + t26 * t50
               + t270);
        dA[6] = t271
            * (t0 * (t116 * t59 + t273) - t276 * t279 - t282 * t284
               + t282 * t309 + t285 * t289 - t293 * t306 + t313);
        dA[7] = t271
            * (-t1 * t50 + t104 * t277 * t287 * t292 * t304 * t319
               + t104 * t278 * t292 * t304 * t324
               + 2 * t277 * t286 * t287 * t319 - t279 * t321 - t284 * t319
               - 2 * t314 - t327);
        dA[8] = t271
            * (t1 * t51 - t279 * (-t21 - t333) - t289 * t337 - t330 + t332
               - t341 - t4 * (-t235 - t328));
        dA[9] = t271
            * (t104 * t277 * t278 * t304 * t350 * t352
               + t104 * t278 * t292 * t304 * t356
               + 2 * t277 * t300 * t350 * t352 - t307 * t347 - t344
               - t350 * t351 - t355);
        dA[10] = t271
            * (t1 * t52 - t2 * (-t26 - t343) - t307 * (-t359 - t41) + t357
               + t358 - t362 + t363 * t365 - t370 + t371 * t372);
        dA[11] = t271
            * (-t1 * t49 + t104 * t278 * t292 * t304 * t381 + t283 * t292 * t376
               - t307 * t374 - 2 * t315 - t365 * t377 - t376 * t382 - t384);
        dA[12] = t48 * (t131 * t388 + t131 * t390 + t175 + t386 * t80);
        dA[13] = t48
            * (-t107
                   * (4 * t145 * t168 + t394 * t395 + t394 * t396
                      + t70 * (t392 + t95) + t77 * (t393 + t97) + t78)
               + t139 + std::pow(t146, 2) * t93 + t146 * t386 + t147 * t388
               + t147 * t390 + 2 * t162 * t397 + t164 * t391 + t399);
        dA[14] = t48 * (t179 * t386 + t179 * t400 + t181 * t387 + t410);
        dA[15] = t48
            * (t0 * (t233 + t413)
               - t107
                   * (-t173 - t196 * t406 - t196 * t418 - t197 * t404
                      - t197 * t407 + t69 * t7 * (t416 + t417)
                      + t7 * t76 * (t144 * t37 + t256))
               + t142 * t35 - t148 - t198 * t386 - t198 * t400 - t198 * t402
               - t198 * t409 - t198 * t415 + t2 * (t220 + t236 + t414) - t200
               - t210 * t397 + t26 * t37 + t4 * (-t17 * t50 + t29 * t52));
        dA[16] = t48
            * (-t107
                   * (2 * t143 * t239 * t7 + 2 * t145 * t238 * t7 - t164 * t395
                      - t164 * t396 + 2 * t2 * t238 * t69 + 2 * t2 * t239 * t76
                      - t70 * (t212 + t423) - t77 * (t214 + t424) - t78)
               + t138 - t226 + t240 * t386 + t240 * t402 + t240 * t409
               + t240 * t415 + t251 - t252 * t397 + t398 + t419 - t421 - t422);
        dA[17] = t48
            * (t0 * (t12 * t144 - t142 * t32)
               - t107
                   * (t259 * t406 + t259 * t418 + t260 * t404 + t260 * t407
                      + t405 + t70 * (t26 + t52 * t9) + t77 * (t429 + t45))
               + t15 * t37 + t153 + t2 * (t257 + t426 + t428) + t261 * t386
               + t261 * t400 + t261 * t402 + t261 * t409 + t261 * t415
               - t269 * t397 + t4 * (t233 + t420) + t401 + t425);
        dA[18] = t271
            * (-t279 * t434 - t282 * t432 + t282 * t441 + t285 * t436
               - t293 * t440 + t3 * t50 + t430 + t443);
        dA[19] = t271
            * (t104 * t277 * t287 * t292 * t319 * t438
               + t104 * t278 * t292 * t324 * t438
               + 2 * t277 * t287 * t319 * t435 - t279 * t445 - t319 * t432
               - t446);
        dA[20] = t271
            * (-t144 * t3 + 2 * t2 * t41 + t278 * t336 * t431 - t279 * t448
               - t336 * t441 - t337 * t436 - t340 * t440
               + t4 * (t50 * t55 + t68) - t447);
        dA[21] = t271
            * (t0 * (t142 * t66 + t449)
               + t104 * t277 * t278 * t350 * t352 * t438
               + t104 * t278 * t292 * t356 * t438
               + 2 * t277 * t350 * t352 * t437 - t3 * t52 - t307 * t451
               - 2 * t335 - t350 * t452 - t357);
        dA[22] = t271
            * (t2 * (t157 * t55 + t454) - t307 * t456 - t361 * t452
               + t363 * t457 - t368 * t458 + t371 * t459 + t460);
        dA[23] = t271
            * (t104 * t278 * t292 * t381 * t438 + t142 * t3 + t292 * t376 * t431
               - t307 * t464 - t377 * t457 - t4 * (-t333 - t463) + t461 - t462
               - t466);
        dA[24] = t48 * (-t131 * t470 + t131 * t471 + t195 - t468 * t80);
        dA[25] = t48 * (-t146 * t468 - t147 * t470 + t147 * t471 + t410);
        dA[26] = t48
            * (-t107
                   * (-t176 * t472 * t69 + 4 * t177 * t403 - t177 * t472 * t76
                      + t70 * (t392 + t94) + t77 * (t393 + t96) + t78)
               + t140 + std::pow(t179, 2) * t93 - t179 * t468 + t179 * t474
               + t180 * t473 - t181 * t469 + t184 * t192 + t399);
        dA[27] = t48
            * (2 * t102 * t104 * t176 * t198 * t83 * t91
               + 2 * t103 * t104 * t177 * t198 * t82 * t91
               - t107
                   * (t194 - t196 * t480 - t197 * t481 + t197 * t483
                      + t467 * t482 + t70 * (t142 * t17 + t463)
                      + t77 * (t232 + t479))
               + t144 * t15 - t179 * t198 * t93 - t185 - t198 * t473
               - t198 * t474 + t2 * t258 - t203 - t210 * t478 + t26 * t29
               + t4 * (t217 + t236 + t477) - t476);
        dA[28] = t48
            * (t102 * t103 * t104 * t190 * t240
               + 2 * t102 * t103 * t240 * t4 * t81 * t91
               - t107
                   * (t238 * t480 - t238 * t486 + t239 * t481 - t239 * t483
                      + t408 + t70 * (t484 + t485) + t77 * (t235 + t29 * t51))
               + t179 * t240 * t81 * t82 * t83 * t91 - t240 * t468 - t240 * t487
               - t252 * t478 - t489);
        dA[29] = t48
            * (t102 * t103 * t104 * t190 * t261
               - t107
                   * (2 * t176 * t4 * t69 + 2 * t177 * t4 * t76
                      + 2 * t259 * t4 * t69 - t259 * t486 + 2 * t260 * t4 * t76
                      - t260 * t483 - t70 * (t211 + t423) - t77 * (t213 + t424)
                      - t78)
               + t179 * t261 * t81 * t82 * t83 * t91 - t261 * t468 - t261 * t487
               - t269 * t478 + t490 - t492);
        dA[30] = t271
            * (t0 * (t144 * t59 + t493)
               + t104 * t277 * t282 * t287 * t292 * t499
               + 2 * t277 * t282 * t287 * t497 - t279 * t494 - t282 * t496
               - t293 * t501 - t329 + 2 * t4 * t45 - t5 * t51);
        dA[31] = t271
            * (t144 * t5 - t2 * (-t256 - t359) - t279 * t503 + t447 + t502
               + t504 * t505 + t508);
        dA[32] = t271
            * (t278 * t336 * t495 - t279 * t510 - t336 * t507 - t337 * t505
               - t340 * t501 + t4 * (t183 * t55 + t509) - t512);
        dA[33] = t271
            * (-t307 * t515 - t350 * t514 + t350 * t520 + t49 * t5 + t513
               + t516 * t517 + t518 * t519 + t521);
        dA[34] = t271
            * (t104 * t277 * t278 * t352 * t361 * t499 - t142 * t5
               + t2 * (t49 * t55 + t61) + 2 * t21 * t4
               + 2 * t277 * t352 * t361 * t498 - t307 * t522 - t361 * t514
               - t368 * t519 - t461);
        dA[35] = t271
            * (t104 * t278 * t292 * t381 * t499 + t292 * t376 * t495
               - t307 * t523 - t376 * t520 - t377 * t517 - t524);
        dA[36] = t271
            * (2 * t104 * t277 * t278 * t304 * t352 * t525
               + 2 * t104 * t277 * t287 * t292 * t304 * t528 - t228
               + t277 * t287 * t304 * t352 * t540 * t81 - t284 * t541
               - t369 * t539
               - t537
                   * (t286 * t527 + t286 * t531 + t299 * (t114 * t9 + t52 * t72)
                      + t300 * t529 + t300 * t532
                      + t302 * (t116 * t29 + t121 * t51) + t533 * t54
                      + t534 * t54 + t535)
               - t54 * t544);
        dA[37] = t271
            * (t0 * (-t155 * t72 + t411 + t546) + t151
               + t2 * (t121 * t157 + t414 + t547) - t200 + t201
               + t4 * (-t157 * t29 + t50 * t72) - t432 * t541 + t438 * t550
               - t438 * t551 + t441 * t552 - t458 * t539 + t465 * t553
               - t537
                   * (t164 * t533 + t164 * t534 + t299 * (t416 + t548)
                      + t302 * (t121 * t155 + t493) + t435 * t527 + t435 * t531
                      + t437 * t529 + t437 * t532));
        dA[38] = t271
            * (2 * t104 * t277 * t278 * t352 * t499 * t525
               + 2 * t104 * t277 * t287 * t292 * t499 * t528 - t186
               + t2 * (t121 * t49 - t183 * t9) - t204
               + t277 * t287 * t352 * t499 * t540 * t81
               + t4 * (t183 * t72 + t477 + t547) - t476 - t496 * t541
               - t499 * t551 - t519 * t539
               - t537
                   * (t192 * t533 + t192 * t534 + t299 * (t189 * t72 + t449)
                      + t302 * (t479 + t545) + t497 * t527 + t497 * t531
                      + t498 * t529 + t498 * t532));
        dA[39] = t48
            * (-t107
                   * (t196 * t69 * t99 + 4 * t197 * t482 + t197 * t76 * t99
                      + t70 * (t560 + t561) + t77 * (t562 + t563) + t78)
               + t129 * t197 * t558 + t196 * t558 * t559
               + std::pow(t198, 2) * t93 + t198 * t210 * t557 + t199 + 2 * t222
               + 2 * t223 - t54 * t554 - t555 - t556);
        dA[40] = t271 * (-t551 * t566 + t552 * t568 + t553 * t569 + t575);
        dA[41] = t271 * (-t551 * t578 + t552 * t579 + t553 * t581 + t585);
        dA[42] = t271
            * (t0 * (-t121 * t4 + t29 * t59) + t279 * t588 - t282 * t539
               + t282 * t592 + t285 * t589 - t293 * t591 + t311 + t511);
        dA[43] = t271
            * (t104 * t277 * t287 * t292 * t319 * t540
               + t104 * t278 * t292 * t324 * t540 + t115 * t55
               + 2 * t277 * t287 * t319 * t528 - t279 * t594 + t29 * t297
               - t296 * t63 - t30 * t55 - t319 * t539 - t593);
        dA[44] = t271
            * (-t1 * t37 + t278 * t336 * t538 - t279 * t595 - t332 - t336 * t592
               - t337 * t589 - t340 * t591 - t596);
        dA[45] = t271
            * (t307 * t600 - t350 * t601 + t350 * t603 + t355 + t516 * t602
               + t518 * t541 - t597);
        dA[46] = t271
            * (-t1 * t17 + t104 * t277 * t278 * t352 * t361 * t540
               + 2 * t277 * t352 * t361 * t525 - t307 * t604 - t358
               - t361 * t601 - t368 * t541 - t605);
        dA[47] = t271
            * (-t10 * t55 + t104 * t278 * t292 * t381 * t540 + t117 * t55
               + t292 * t376 * t538 + t296 * t9 - t297 * t56 - t307 * t607
               - t376 * t603 - t377 * t602 - t606);
        dA[48] = t271
            * (2 * t104 * t277 * t278 * t304 * t352 * t565
               + 2 * t104 * t277 * t287 * t292 * t304 * t564 - t110
               - t164 * t544 + t229 * t4 - t234 - t237
               + t277 * t287 * t304 * t352 * t566 * t81 - t284 * t614
               - t369 * t572
               - t537
                   * (t286 * t609 + t286 * t612 + t299 * (t114 * t58 + t61)
                      + t300 * t608 + t300 * t613 + t302 * (t124 + t241)
                      + t54 * t610 + t54 * t611)
               - t574);
        dA[49] = t271
            * (2 * t104 * t277 * t278 * t352 * t438 * t565
               + 2 * t104 * t277 * t287 * t292 * t438 * t564 - t208 - t223
               - t227 + t277 * t287 * t352 * t438 * t566 * t81 - t421
               - t432 * t614 - t438 * t616 - t458 * t572 - t491
               - t537
                   * (t164 * t610 + t164 * t611
                      + t299 * (t142 * t58 + t157 * t17)
                      + t302 * (t155 * t32 + t50 * t63) + t435 * t609
                      + t435 * t612 + t437 * t608 + t437 * t613 + t535));
        dA[50] = t271
            * (2 * t104 * t277 * t278 * t352 * t499 * t565
               + 2 * t104 * t277 * t287 * t292 * t499 * t564
               + t277 * t287 * t352 * t499 * t566 * t81 - t489 - t496 * t614
               - t499 * t616 - t519 * t572
               - t537
                   * (t192 * t610 + t192 * t611 + t299 * (t484 + t617)
                      + t302 * (t183 * t63 + t326) + t497 * t609 + t497 * t612
                      + t498 * t608 + t498 * t613));
        dA[51] =
            t271 * (-t540 * t616 + 2 * t564 * t592 + 2 * t565 * t603 + t575);
        dA[52] = t48
            * (2 * t102 * t104 * t238 * t240 * t83 * t91
               + 2 * t103 * t104 * t239 * t240 * t82 * t91
               - t107
                   * (t238 * t239 * t621 - t238 * t394 * t69 - t239 * t394 * t76
                      + t70 * (t561 + t619) + t77 * (t563 + t620) + t78)
               - t164 * t618 + std::pow(t240, 2) * t81 * t82 * t83 * t91
               - t240 * t252 * t557 + 2 * t29 * t61 - t419 - t556 - t622);
        dA[53] = t271
            * (2 * t104 * t277 * t278 * t352 * t565 * t578
               + 2 * t104 * t277 * t287 * t292 * t564 * t578 - t578 * t616
               - t626);
        dA[54] = t271
            * (t104 * t277 * t282 * t287 * t292 * t566
               + 2 * t277 * t282 * t287 * t564 - t279 * t628 - t282 * t572
               - t293 * t629 - t430 - t627 - t631);
        dA[55] = t271
            * (t2 * t35 + t279 * t633 + t319 * t568 - t319 * t572 + t323
               + t504 * t634 + t506 * t629);
        dA[56] = t271
            * (t124 * t59 + t278 * t336 * t571 - t279 * t636 + t295 * t32
               - t297 * t65 - t336 * t568 - t337 * t634 - t340 * t629
               - t39 * t59 - t635);
        dA[57] = t271
            * (t104 * t277 * t278 * t350 * t352 * t566
               + t104 * t278 * t292 * t356 * t566 + t156 * t59 + t17 * t297
               - t23 * t59 + 2 * t277 * t350 * t352 * t565 - t307 * t639 - t637
               - t638 - t641);
        dA[58] = t271
            * (t104 * t277 * t278 * t352 * t361 * t566 + t277 * t292 * t643
               + 2 * t277 * t352 * t361 * t565 - t361 * t640 - t368 * t614
               - t644);
        dA[59] = t271
            * (t104 * t278 * t292 * t381 * t566 + t292 * t376 * t571
               - t307 * t646 - t364 * t377 * t565 - t376 * t569 + t462 - t645
               - t648);
        dA[60] = t271
            * (t0 * (t114 * t65 + t255 + t649) + t112 - t192 * t544
               + t2 * (-t114 * t37 + t51 * t56) + t264 + t270 - t284 * t624
               + t309 * t657 - t369 * t584 + t382 * t658
               - t537
                   * (t286 * t650 + t286 * t654 + t299 * (t120 + t262)
                      + t300 * t651 + t300 * t655 + t302 * (t116 * t65 + t68)
                      + t54 * t652 + t54 * t653)
               + t542 * t656);
        dA[61] = t271
            * (t0 * (-t12 * t155 + t142 * t65) + t154
               + t2 * (t155 * t56 + t428 + t649)
               + t4 * (-t157 * t65 + t231 + t546) - t432 * t624
               + t438 * t656 * t81 - t438 * t659 + t441 * t657 - t458 * t584
               + t465 * t658
               - t537
                   * (t164 * t652 + t164 * t653 + t299 * (t157 * t56 + t75)
                      + t302 * (t160 + t429) + t435 * t650 + t435 * t654
                      + t437 * t651 + t437 * t655)
               + t625);
        dA[62] = t271
            * (2 * t104 * t277 * t278 * t352 * t499 * t576
               + 2 * t104 * t277 * t287 * t292 * t499 * t577
               + t277 * t287 * t352 * t499 * t578 * t81 - t492 - t496 * t624
               - t499 * t659 - t519 * t584
               - t537
                   * (t192 * t652 + t192 * t653
                      + t299 * (t12 * t189 + t49 * t56)
                      + t302 * (t144 * t65 + t183 * t37) + t497 * t650
                      + t497 * t654 + t498 * t651 + t498 * t655 + t535));
        dA[63] = t271 * (-t540 * t659 + t585 + t592 * t657 + t603 * t658);
        dA[64] = t271
            * (2 * t104 * t277 * t278 * t352 * t566 * t576
               + 2 * t104 * t277 * t287 * t292 * t566 * t577 - t566 * t659
               - t626);
        dA[65] = t48
            * (2 * t102 * t104 * t259 * t261 * t83 * t91
               + 2 * t103 * t104 * t260 * t261 * t82 * t91
               - t107
                   * (t259 * t260 * t621 - t259 * t472 * t69 - t260 * t472 * t76
                      + t70 * (t560 + t619) + t77 * (t562 + t620) + t78)
               - t192 * t623 + std::pow(t261, 2) * t81 * t82 * t83 * t91
               - t261 * t269 * t557 - t490 - t555 - t622 + 2 * t68 * t9);
        dA[66] = t271
            * (t160 * t66 + 2 * t277 * t282 * t287 * t577 - t279 * t660
               - t43 * t66 - t663);
        dA[67] = t271
            * (t104 * t277 * t287 * t292 * t319 * t578
               + t104 * t278 * t292 * t324 * t578
               + 2 * t277 * t287 * t319 * t577 - t279 * t665 - t319 * t584
               - t502 - t664 - t667);
        dA[68] = t271
            * (t279 * t668 - t288 * t337 * t577 - t336 * t579 + t336 * t584
               - t340 * t662 + t4 * (-t2 * t65 + t37 * t55) + t512);
        dA[69] = t271
            * (t104 * t277 * t278 * t350 * t352 * t578
               + t104 * t278 * t292 * t356 * t578
               + 2 * t277 * t350 * t352 * t576 - t307 * t671 - t350 * t670
               - t513 - t669 - t673);
        dA[70] = t271
            * (t104 * t277 * t278 * t352 * t361 * t578 + t12 * t295 + t120 * t66
               - t19 * t66 + 2 * t277 * t352 * t361 * t576 - t296 * t58
               - t307 * t675 - t361 * t670 - t368 * t624 - t674);
        dA[71] = t271
            * (t104 * t278 * t292 * t381 * t578 + t277 * t292 * t676
               + t292 * t376 * t583 - t364 * t377 * t576 - t376 * t581 - t677);
        dA[72] = t48
            * (t0 * (t273 + t328) - t107 * (t276 * t70 + t683 * t71 + t684)
               + t126 * t680 + t131 * t685 + t131 * t686 + t312 + t678
               - t682 * t80);
        dA[73] = t48
            * (-t107 * (t145 * t683 + t434 * t70 + t687) - t146 * t682
               + t146 * t688 + t147 * t686 + t162 * t680 - t2 * t464 + t375
               + t443);
        dA[74] = t48
            * (t0 * (-t41 - t450)
               - t107 * (-t177 * t683 + t282 * t480 + t494 * t70) - t179 * t682
               + t179 * t688 + t179 * t689 + t190 * t680 + t330 + t4 * t45
               + t4 * t522);
        dA[75] = t48
            * (-t0 * t45 + t107 * (t197 * t683 + t588 * t70 + t684)
               + t198 * t682 - t198 * t688 - t198 * t689 - t210 * t680 + t311
               + t690);
        dA[76] = t48
            * (2 * t102 * t104 * t240 * t282 * t83 * t91
               + t102 * t104 * t240 * t679 * t83 * t91
               - t107 * (t239 * t683 + t628 * t70 - t687) - t2 * t646
               - t240 * t682 - t252 * t680 - t375 - t631);
        dA[77] = t271
            * (2 * t104 * t277 * t287 * t292 * t578 * t692
               - t537 * (t299 * (t19 + t449) + t650 * t692 + t654 * t692)
               + t66 * (-t326 - t43) - t663);
        dA[78] = t695
            * (t279 * std::pow(t282, 2) - t285 * t291 + t285 * t694 + t693);
        dA[79] = t695 * (t504 * t694 + t699);
        dA[80] = t695 * (t291 * t336 - t337 * t694 - t700 - t701 - t702);
        dA[81] = t706;
        dA[82] = t709;
        dA[83] = t712;
        dA[84] = t48
            * (t0 * (t0 * t116 + t545) + t102 * t103 * t104 * t324 * t80
               + 2 * t102 * t104 * t319 * t80 * t83 * t91
               - t107 * (t321 * t70 + t71 * t714 + t715) - t126 * t713
               - t131 * t716 - t314 - t327);
        dA[85] = t48
            * (t102 * t103 * t104 * t146 * t324
               + 2 * t102 * t104 * t146 * t319 * t83 * t91
               - t107 * (t145 * t714 + t445 * t70 + t717) - t146 * t718
               - t162 * t713 - t446);
        dA[86] = t271
            * (-t455 * t55 + 2 * t507 * t720 + t508
               - t537
                   * (t192 * t294 * t720 + t299 * (t49 * t66 + t75)
                      + t498 * t526 * t720)
               + t59 * (t0 * t183 + t124) + t66 * (t144 * t66 + t68) - t719);
        dA[87] = t48
            * (t0 * t35 + t0 * t607
               - t107 * (-t197 * t714 + t594 * t69 * t7 - t715) + t198 * t718
               - t198 * t722 - t198 * t723 + t210 * t713 + t627 + t630);
        dA[88] = t48
            * (-t107 * (2 * t239 * t319 * t7 - t633 * t70 - t717) - t240 * t718
               + t240 * t722 + t240 * t723 + t252 * t713 + t323 + t678);
        dA[89] = t48
            * (-t107 * (t260 * t714 - t319 * t480 + t665 * t70) - t261 * t718
               + t261 * t722 + t261 * t723 + t269 * t713 + t349
               + t4 * (t64 - 2 * t67) - t664 - t666);
        dA[90] = t695 * (t285 * t724 + t699);
        dA[91] = t695
            * (t279 * std::pow(t319, 2) + t324 * t504 + t504 * t724 + t725);
        dA[92] = t695 * (-t337 * t724 - t727);
        dA[93] = t729;
        dA[94] = t731;
        dA[95] = t733;
        dA[96] = t271
            * (2 * t104 * t277 * t287 * t292 * t304 * t738 - t341 - t345 * t59
               - t537 * (t294 * t54 * t738 + t299 * t736 + t300 * t526 * t738)
               + t55 * t735 + t66 * (t116 * t2 + t160) - t734);
        dA[97] = t48
            * (-t107 * (-t145 * t745 + t448 * t69 * t7 - t744)
               - t132 * t147 * t746 - t146 * t743 - t162 * t740 + t2 * t451
               + t4 * (-t232 - t373) - t447 - t739);
        dA[98] = t48
            * (-t107 * (t177 * t745 + t510 * t70 + t747) - t179 * t743
               - t179 * t748 - t190 * t740 - t322 + t4 * (t359 + t509) - t690);
        dA[99] = t48
            * (-t0 * t604 + t102 * t103 * t104 * t198 * t741
               + 3 * t102 * t104 * t198 * t336 * t83 * t91 + t103 * t210 * t336
               - t107 * (t171 * t336 + t197 * t745 + t595 * t70) - t331 - t596);
        dA[100] = t48
            * (-t107 * (-t239 * t745 + t636 * t70 + t744) + t2 * t639
               - t240 * t743 - t240 * t748 + t252 * t740 + t667 + t739);
        dA[101] = t48
            * (t103 * t269 * t336 - t107 * (-t260 * t745 - t668 * t70 - t747)
               - t261 * t743 - t261 * t748 + t322 - t4 * t41 - t678);
        dA[102] = t695 * (t285 * t750 + t291 * t336 - t700 - t702 + t749);
        dA[103] = t695 * (2 * t277 * t278 * t319 * t738 - t727);
        dA[104] = t695
            * (t277 * t278 * std::pow(t336, 2) + 2 * t336 * t339 - t337 * t750
               - t751);
        dA[105] = t753;
        dA[106] = t754;
        dA[107] = t755;
        dA[108] = t48
            * (-t107 * (t347 * t77 + t516 * t98 + t756) + t126 * t704
               + t127 * t356 - t344 + t757 * t80 + t758 * t80 + t759);
        dA[109] = t48
            * (t0 * (-t21 - t453)
               - t107 * (t168 * t516 + t350 * t407 + t451 * t77) + t146 * t757
               + t146 * t758 + t162 * t704 - t2 * t26 + t2 * t448 + t356 * t397
               - t357);
        dA[110] = t48
            * (-t107 * (-t403 * t516 + t515 * t77 + t760) + t179 * t757
               + t179 * t758 + t190 * t704 + t317 + t356 * t478 - t4 * t503
               + t521);
        dA[111] = t48
            * (t102 * t103 * t104 * t91 * (t482 * t516 + t600 * t77 + t756)
               - t198 * t757 - t198 * t758 - t198 * t761 - t210 * t704 - t597
               - t759);
        dA[112] = t271
            * (t17 * t297 + t350 * t569 + t518 * t614
               - t537 * (t302 * (t39 + t493) + t608 * t762 + t613 * t762)
               + 2 * t569 * t762 + t59 * (-t23 - t75) - t637 - t638 - t641);
        dA[113] = t48
            * (-t107 * (t516 * t763 + t671 * t77 - t760) + t261 * t757
               + t261 * t758 + t261 * t761 - t269 * t704 + t316
               + t4 * (t73 - 2 * t74) - t669 - t672);
        dA[114] = t706;
        dA[115] = t729;
        dA[116] = t753;
        dA[117] = t765
            * (t307 * std::pow(t350, 2) + t356 * t516 + t516 * t764 + t693);
        dA[118] = t765 * (t363 * t764 + t767);
        dA[119] = t765 * (t350 * t381 - t377 * t764 - t701 - t768 - t769);
        dA[120] = t271
            * (t104 * t277 * t278 * t304 * t352 * t361
               + 2 * t104 * t277 * t278 * t304 * t352 * t770 - t275 * t66 - t362
               - t370
               - t537 * (t286 * t526 * t770 + t301 * t54 * t770 + t302 * t735)
               + t55 * t736 + t59 * (t114 * t4 + t617) - t737);
        dA[121] = t48
            * (t0 * t26 - t107 * (t168 * t363 + t456 * t77 + t771) + t146 * t772
               - t146 * t773 + t15 * t4 - t162 * t708 + t2 * (t333 + t454)
               - t397 * t707);
        dA[122] = t48
            * (2 * t103 * t104 * t179 * t361 * t82 * t91
               - t107 * (-t363 * t403 + t522 * t77 + t775) - t179 * t773
               - t190 * t708 + t2 * (-t320 - t417) + t4 * t494 - t461
               - t478 * t707 - t774);
        dA[123] = t48
            * (-t0 * t595 + t102 * t103 * t104 * t198 * t707
               + t102 * t210 * t360 + t103 * t104 * t198 * t360 * t82 * t91
               - t107 * (-t172 * t361 - t363 * t482 + t604 * t7 * t76)
               - t198 * t772 - t334 - t605);
        dA[124] = t48
            * (t102 * t252 * t360 + 2 * t103 * t104 * t240 * t361 * t82 * t91
               - t107 * (2 * t238 * t361 * t7 - t643 * t77 - t771) - t240 * t773
               - t240 * t776 - t644);
        dA[125] = t48
            * (-t107 * (t363 * t763 + t675 * t77 - t775) + t261 * t772
               - t261 * t773 - t261 * t776 + t269 * t708 + t4 * t660 + t645
               + t647 + t774);
        dA[126] = t709;
        dA[127] = t731;
        dA[128] = t754;
        dA[129] = t765 * (t516 * t777 + t767);
        dA[130] = t765
            * (t307 * std::pow(t361, 2) - t363 * t367 + t363 * t777 + t725);
        dA[131] = t765 * (-t377 * t777 + t778);
        dA[132] = t48
            * (t0 * (t0 * t114 + t548) + t102 * t103 * t104 * t381 * t80
               - t107 * (t374 * t7 * t76 - t377 * t98 - t779) - t126 * t732
               - t315 - t384 - t780 * t80);
        dA[133] = t271
            * (t278 * t381 * t458 + t376 * t452 - t444 * t55 + 2 * t465 * t781
               - t466
               - t537
                   * (t164 * t301 * t781 + t302 * (t326 + t50 * t59)
                      + t435 * t526 * t781)
               + t59 * (t142 * t59 + t61) + t66 * (t0 * t157 + t120) - t691);
        dA[134] = t48
            * (t102 * t103 * t104 * t179 * t381
               - t107 * (t377 * t403 + t523 * t77 + t782) - t179 * t780
               - t190 * t732 - t524);
        dA[135] = t48
            * (t0 * t15 + t0 * t594 - t106 * t198 * t381
               - t107 * (t377 * t482 + t607 * t77 + t779) + t198 * t780
               + t210 * t732 + t673);
        dA[136] = t48
            * (t102 * t103 * t104 * t240 * t381 + t102 * t252 * t376
               - t107 * (-t238 * t377 * t7 + t376 * t407 + t646 * t77)
               + t2 * t61 + t2 * (t57 - 2 * t60) - t240 * t780 - t648);
        dA[137] = t48
            * (t102 * t103 * t104 * t261 * t381 + t102 * t269 * t376
               - t107 * (-t377 * t763 - t676 * t77 - t782) - t261 * t780
               - t677);
        dA[138] = t712;
        dA[139] = t733;
        dA[140] = t755;
        dA[141] = t765 * (t350 * t381 + t516 * t783 + t749 - t768 - t769);
        dA[142] = t765 * (t363 * t783 + t778);
        dA[143] = t765
            * (t277 * t292 * std::pow(t376, 2) - t377 * t381 - t377 * t783
               - t751);
    }

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
        double dA[9])
    {
        const auto t0 = -t2_y;
        const auto t1 = t0 + t1_y;
        const auto t2 = t0_x - t1_x;
        const auto t3 = t0 + t0_y;
        const auto t4 = -t2_x;
        const auto t5 = t0_x + t4;
        const auto t6 = t0_y - t1_y;
        const auto t7 = t2 * t3 - t5 * t6;
        const auto t8 = -t2_z;
        const auto t9 = t1_z + t8;
        const auto t10 = t0_z + t8;
        const auto t11 = t0_z - t1_z;
        const auto t12 = t10 * t2 - t11 * t5;
        const auto t13 =
            1.0 / (std::pow(t11, 2) + std::pow(t2, 2) + std::pow(t6, 2));
        const auto t14 = t10 * t6 - t11 * t3;
        const auto t15 = std::pow(t14, 2) + std::pow(t7, 2);
        const auto t16 = std::pow(t12, 2) + t15;
        const auto t17 = t10 * t2 - t11 * t5;
        const auto t18 = t15 + std::pow(t17, 2);
        const auto t19 = t13 * t18;
        const auto t20 = std::sqrt(t19) / t18;
        const auto t21 = t1_x + t4;
        const auto t22 = t13 * t16;
        dA[0] = -t20 * (-t1 * t7 - t12 * t9 + t13 * t16 * t2);
        dA[1] = -t20 * (-t14 * t9 + t21 * t7 + t22 * t6);
        dA[2] = -t20 * (t1 * t14 + t11 * t22 + t12 * t21);
        dA[3] = t20 * (-t10 * t12 + t13 * t18 * t2 - t3 * t7);
        dA[4] = t20 * (-t10 * t14 + t19 * t6 + t5 * t7);
        dA[5] = t20 * (t11 * t19 + t14 * t3 + t17 * t5);
        dA[6] = t20 * (t11 * t17 + t6 * t7);
        dA[7] = -t20 * (-t11 * t14 + t2 * t7);
        dA[8] = -t20 * (t12 * t2 + t14 * t6);
    }

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
        double dA[81])
    {
        const auto t0 = -t2_y;
        const auto t1 = t0 + t1_y;
        const auto t2 = std::pow(t1, 2);
        const auto t3 = t0_x - t1_x;
        const auto t4 = -t2_z;
        const auto t5 = t0_z + t4;
        const auto t6 = t3 * t5;
        const auto t7 = -t2_x;
        const auto t8 = t0_x + t7;
        const auto t9 = t0_z - t1_z;
        const auto t10 = t8 * t9;
        const auto t11 = -t10;
        const auto t12 = t11 + t6;
        const auto t13 = t0 + t0_y;
        const auto t14 = t13 * t3;
        const auto t15 = t0_y - t1_y;
        const auto t16 = t15 * t8;
        const auto t17 = t14 - t16;
        const auto t18 = t15 * t5;
        const auto t19 = t13 * t9;
        const auto t20 = -t19;
        const auto t21 = t18 + t20;
        const auto t22 = std::pow(t17, 2) + std::pow(t21, 2);
        const auto t23 = std::pow(t12, 2) + t22;
        const auto t24 = 1.0 / t23;
        const auto t25 = std::pow(t9, 2);
        const auto t26 = std::pow(t3, 2);
        const auto t27 = std::pow(t15, 2);
        const auto t28 = t26 + t27;
        const auto t29 = t25 + t28;
        const auto t30 = 1.0 / t29;
        const auto t31 = -t3;
        const auto t32 = -t5;
        const auto t33 = -t8;
        const auto t34 = -t9;
        const auto t35 = t31 * t32 - t33 * t34;
        const auto t36 = t22 + std::pow(t35, 2);
        const auto t37 = -t3 * t30 * t36;
        const auto t38 = t1 * t17;
        const auto t39 = t1_z + t4;
        const auto t40 = t35 * t39 + t38;
        const auto t41 = t37 + t40;
        const auto t42 = -t41;
        const auto t43 = t12 * t39 + t38;
        const auto t44 = t3 * t30;
        const auto t45 = 2 * t44;
        const auto t46 = t24 * t42;
        const auto t47 = 2 * t40;
        const auto t48 = 4 / std::pow(t29, 2);
        const auto t49 = t23 * t48;
        const auto t50 = t26 * t49;
        const auto t51 = t23 * t30;
        const auto t52 = -t51;
        const auto t53 = std::pow(t39, 2) + t52;
        const auto t54 = t24 * std::sqrt(t51);
        const auto t55 = t30 * t36;
        const auto t56 = t15 * t55;
        const auto t57 = t1_x + t7;
        const auto t58 = t17 * t57 - t21 * t39;
        const auto t59 = t56 + t58;
        const auto t60 = t30 * t47;
        const auto t61 = t15 * t60;
        const auto t62 = -t58;
        const auto t63 = t15 * t3;
        const auto t64 = t36 * t48;
        const auto t65 = t63 * t64;
        const auto t66 = -t65;
        const auto t67 = t45 * t62 + t66;
        const auto t68 = t1 * t57 - t24 * t42 * t59 + t61 + t67;
        const auto t69 = t39 * t57;
        const auto t70 = t55 * t9;
        const auto t71 = t1 * t21 + t35 * t57;
        const auto t72 = t70 + t71;
        const auto t73 = t45 * t72;
        const auto t74 = t45 * t71;
        const auto t75 = t24 * t72;
        const auto t76 = t60 * t9;
        const auto t77 = t3 * t9;
        const auto t78 = t64 * t77;
        const auto t79 = -t76 + t78;
        const auto t80 = t13 * t17 + t35 * t5;
        const auto t81 = t23 * t3 * t30 - t80;
        const auto t82 = t24 * t81;
        const auto t83 = t39 * t5;
        const auto t84 = t52 + t83;
        const auto t85 = t1 * t13;
        const auto t86 = -2 * t3 * t30 * t80 + t85;
        const auto t87 = t17 * t8 - t21 * t5;
        const auto t88 = t15 * t51 + t87;
        const auto t89 = t24 * t88;
        const auto t90 = t45 * t87;
        const auto t91 = -t14;
        const auto t92 = t16 + t91;
        const auto t93 = -t1 * t8 - t61 + t65 + t90 + t92;
        const auto t94 = t13 * t21;
        const auto t95 = t12 * t8 + t94;
        const auto t96 = t51 * t9 + t95;
        const auto t97 = t24 * t96;
        const auto t98 = t35 * t8 + t94;
        const auto t99 = t45 * t98;
        const auto t100 = -t6;
        const auto t101 = t10 + t100;
        const auto t102 = t101 - t39 * t8 + t79 + t99;
        const auto t103 = t15 * t17;
        const auto t104 = t103 + t12 * t9;
        const auto t105 = t24 * t47;
        const auto t106 = t1 * t15;
        const auto t107 = t39 * t9;
        const auto t108 = -t104 * t46 + t106 + t107;
        const auto t109 = t1 * t3;
        const auto t110 = t17 * t3 - t21 * t9;
        const auto t111 = t3 * t39;
        const auto t112 = t15 * t21 + t3 * t35;
        const auto t113 = t112 * t46;
        const auto t114 = t15 * t30;
        const auto t115 = 2 * t114;
        const auto t116 = std::pow(t57, 2);
        const auto t117 = t24 * t59;
        const auto t118 = 2 * t62;
        const auto t119 = t27 * t49;
        const auto t120 = 4 * t114;
        const auto t121 = t1 * t39;
        const auto t122 = t30 * t9;
        const auto t123 = t118 * t122;
        const auto t124 = t15 * t9;
        const auto t125 = t124 * t64;
        const auto t126 = -t123 + t125;
        const auto t127 = t115 * t71;
        const auto t128 = -t115 * t72 + t127;
        const auto t129 = t115 * t80;
        const auto t130 = t129 + t13 * t57 + t67 + t92;
        const auto t131 = t57 * t8;
        const auto t132 = t115 * t87 + t131;
        const auto t133 = t115 * t98;
        const auto t134 = -t18 + t19;
        const auto t135 = t126 - t13 * t39 + t133 + t134;
        const auto t136 = t118 * t24;
        const auto t137 = t104 * t117 + t15 * t57 + t92;
        const auto t138 = t3 * t57;
        const auto t139 = t107 + t110 * t117 + t138;
        const auto t140 = t15 * t39;
        const auto t141 = t112 * t117;
        const auto t142 = 2 * t122;
        const auto t143 = 2 * t71;
        const auto t144 = -t78;
        const auto t145 = t144 - t74;
        const auto t146 = -t142 * t72;
        const auto t147 = 4 * t122;
        const auto t148 = t25 * t49 + t52;
        const auto t149 = t142 * t80;
        const auto t150 = t101 + t145 + t149 + t5 * t57;
        const auto t151 = t125 + t142 * t87;
        const auto t152 = -t1 * t5 + t151 + t21;
        const auto t153 = t131 + t142 * t71 + t85;
        const auto t154 = t101 + t104 * t75 + t57 * t9;
        const auto t155 = t1 * t9;
        const auto t156 = t143 * t24;
        const auto t157 = t112 * t75;
        const auto t158 = t106 + t138;
        const auto t159 = t26 * t64;
        const auto t160 = 1.0 / t36;
        const auto t161 = -t37 - t80;
        const auto t162 = -t13;
        const auto t163 = -t15;
        const auto t164 = t162 * t31 - t163 * t33;
        const auto t165 = -t35;
        const auto t166 = 2 * t162 * t164 + 2 * t165 * t5;
        const auto t167 = t160 * t166;
        const auto t168 = -t55;
        const auto t169 = t168 + t83;
        const auto t170 = t160 * std::sqrt(t55);
        const auto t171 = -t59;
        const auto t172 = t160 * t161;
        const auto t173 = std::pow(t13, 2);
        const auto t174 = t168 + std::pow(t5, 2);
        const auto t175 = t56 + t87;
        const auto t176 = t160 * t175;
        const auto t177 = t129 + t13 * t8 - t160 * t161 * t175 + t66 - t90;
        const auto t178 = t70 + t98;
        const auto t179 = t160 * t178;
        const auto t180 = t144 + t149 - t160 * t161 * t178 + t5 * t8 - t99;
        const auto t181 = t13 * t15;
        const auto t182 = t5 * t9;
        const auto t183 = -t104 * t24 * t81 + t181 + t182;
        const auto t184 = 2 * t24;
        const auto t185 = t184 * t80;
        const auto t186 = t110 * t82 - 2 * t13 * t3 + t16;
        const auto t187 = t112 * t82;
        const auto t188 = -t162 * t34 + t163 * t32;
        const auto t189 = t164 * t8 + t188 * t32;
        const auto t190 = 2 * t189;
        const auto t191 = t160 * t190;
        const auto t192 = t27 * t64;
        const auto t193 = std::pow(t8, 2);
        const auto t194 = -t13 * t5 + t133 + t151 + t176 * t178;
        const auto t195 = t184 * t87;
        const auto t196 = t3 * t8;
        const auto t197 = t110 * t89 + t182 + t196;
        const auto t198 = t112 * t89;
        const auto t199 = 2 * t18 + t20;
        const auto t200 = t13 * t188 + t165 * t33;
        const auto t201 = 2 * t200;
        const auto t202 = t160 * t201;
        const auto t203 = t168 + t25 * t64;
        const auto t204 = t110 * t97 - 2 * t13 * t9 + t18;
        const auto t205 = t112 * t97;
        const auto t206 = t181 + t196;
        const auto t207 = t103 + t35 * t9;
        const auto t208 = 2 * t207;
        const auto t209 = -t110;
        const auto t210 = t160 * t209;
        const auto t211 = t15 * t164 + t165 * t34;
        const auto t212 = -t160 * t207 * t209 + t63;
        const auto t213 = t160 * t207;
        const auto t214 = t112 * t213 + t77;
        const auto t215 = 2 * t209;
        const auto t216 = t164 * t31 + t188 * t9;
        const auto t217 = t112 * t210 + t124;
        const auto t218 = -t112 * t45;
        const auto t219 = -t112 * t115;
        const auto t220 = t112 * t142;
        const auto t221 = 2 * t163 * t188 + 2 * t165 * t3;
        dA[0] = t54
            * (t2 + t24 * std::pow(t42, 2) - t42 * t45 - 4 * t43 * t44
               + t46 * t47 + t50 + t53);
        dA[1] = t54 * (2 * t24 * t40 * t59 - t45 * t59 - t68);
        dA[2] = t54 * (t46 * t72 + t47 * t75 - t69 - t73 + t74 + t79);
        dA[3] = t54
            * (2 * t3 * t30 * t43 + 2 * t3 * t30 * t81 - t46 * t81 - t47 * t82
               - t50 - t84 - t86);
        dA[4] = t54 * (2 * t3 * t30 * t88 - t46 * t88 - t47 * t89 - t93);
        dA[5] = t54 * (-t102 + 2 * t3 * t30 * t96 - t46 * t96 - t47 * t97);
        dA[6] = t54 * (-t104 * t105 + t108);
        dA[7] = t54 * (t105 * t110 - t109 + t110 * t46 + t92);
        dA[8] = t54 * (t101 + t105 * t112 - t111 + t113);
        dA[9] = t54 * (-t115 * t42 + 2 * t24 * t42 * t62 - t68);
        dA[10] = t54
            * (-t115 * t59 + t116 + t117 * t118 + t119 + t120 * t58
               + t24 * std::pow(t59, 2) + t53);
        dA[11] = t54 * (t117 * t72 + t118 * t75 - t121 + t126 + t128);
        dA[12] = t54 * (t115 * t81 - t117 * t81 - t118 * t82 + t130);
        dA[13] = t54
            * (-t115 * t58 - t117 * t88 - t118 * t89 - t119 - t132
               + 2 * t15 * t30 * t88 - t84);
        dA[14] = t54 * (-t117 * t96 - t118 * t97 - t135 + 2 * t15 * t30 * t96);
        dA[15] = -t54 * (t104 * t136 + t137);
        dA[16] = t54 * (t110 * t136 + t139);
        dA[17] = t54 * (t112 * t136 + t134 - t140 + t141);
        dA[18] = t54
            * (-t142 * t42 - t143 * t46 - t145 + t24 * t42 * t72 - t69 - t76);
        dA[19] = t54
            * (-t117 * t143 - t121 - t123 + t125 + t127 - t142 * t59
               + t24 * t59 * t72);
        dA[20] = t54
            * (t116 - t143 * t75 + t146 + t147 * t71 + t148 + t2
               + t24 * std::pow(t72, 2));
        dA[21] = t54 * (t142 * t81 + t143 * t82 + t150 - t75 * t81);
        dA[22] = t54
            * (-t127 - t152 + 2 * t24 * t71 * t88 + 2 * t30 * t88 * t9
               - t75 * t88);
        dA[23] = t54
            * (-t142 * t95 - t148 - t153 + 2 * t24 * t71 * t96
               + 2 * t30 * t9 * t96 - t75 * t96);
        dA[24] = t54 * (2 * t104 * t24 * t71 - t154);
        dA[25] = t54 * (-t110 * t156 + t110 * t24 * t72 - t134 - t155);
        dA[26] = t54 * (-t112 * t156 + t157 + t158);
        dA[27] = t170
            * (-t159 + t160 * t161 * t41 - t167 * t41 - t169
               + 2 * t3 * t30 * t40 - t41 * t45 - t86);
        dA[28] = t170 * (t130 - t167 * t171 + t171 * t172 - t171 * t45);
        dA[29] = t170 * (t150 + t167 * t72 - t172 * t72 + t73);
        dA[30] = t170
            * (t159 + t160 * std::pow(t161, 2) - t161 * t45 - t166 * t172 + t173
               + t174 - 4 * t44 * t80);
        dA[31] = t170 * (-t166 * t176 - t175 * t45 - t177);
        dA[32] = t170 * (-t166 * t179 - t178 * t45 - t180);
        dA[33] = t54 * (2 * t104 * t24 * t80 - t183);
        dA[34] = t54 * (-t110 * t185 - t186);
        dA[35] = t54 * (-t10 - t112 * t185 - t187 + 2 * t3 * t5);
        dA[36] = t170 * (-t115 * t41 + t160 * t175 * t41 - t191 * t41 - t93);
        dA[37] = t170
            * (-t115 * t171 - t132 + 2 * t15 * t30 * t62 + t160 * t171 * t175
               - t169 - t171 * t191 - t192);
        dA[38] = t170 * (-t128 - t152 + 2 * t160 * t189 * t72 - t176 * t72);
        dA[39] = t170 * (-t115 * t161 - t172 * t190 - t177);
        dA[40] = t170
            * (-t115 * t175 + t120 * t87 + t160 * std::pow(t175, 2) + t174
               - t176 * t190 + t192 + t193);
        dA[41] = t170 * (-t115 * t178 - t179 * t190 + t194);
        dA[42] = t54 * (-t104 * t195 + t104 * t24 * t88 - t14 + 2 * t15 * t8);
        dA[43] = t54 * (2 * t110 * t24 * t87 - t197);
        dA[44] = t54 * (t112 * t195 - t198 + t199);
        dA[45] = t170 * (-t102 - t142 * t41 + t160 * t178 * t41 - t202 * t41);
        dA[46] =
            t170 * (-t135 - t142 * t171 + t160 * t171 * t178 - t171 * t202);
        dA[47] = t170
            * (-t142 * t98 - t146 - t153 + 2 * t160 * t200 * t72 - t179 * t72
               - t203);
        dA[48] = t170 * (-t142 * t161 - t172 * t201 - t180);
        dA[49] = t170 * (-t142 * t175 - t176 * t201 + t194);
        dA[50] = t170
            * (-t142 * t178 + t147 * t98 + t160 * std::pow(t178, 2) + t173
               - t179 * t201 + t193 + t203);
        dA[51] =
            t54 * (-t104 * t184 * t98 + t104 * t24 * t96 - t6 + 2 * t8 * t9);
        dA[52] = t54 * (2 * t110 * t24 * t98 - t204);
        dA[53] = t54 * (2 * t112 * t24 * t98 - t205 - t206);
        dA[54] = t54 * (t108 - t207 * t45 + t208 * t46);
        dA[55] = t54 * (-t115 * t207 - t137 + 2 * t207 * t24 * t59);
        dA[56] = t54 * (-t142 * t207 - t154 + 2 * t207 * t24 * t72);
        dA[57] = t54 * (2 * t104 * t3 * t30 - t183 - t208 * t82);
        dA[58] = t54 * (t104 * t115 + t104 * t89 + 2 * t16 - t208 * t89 + t91);
        dA[59] = t54 * (2 * t10 + t100 + t104 * t142 + t104 * t97 - t208 * t97);
        dA[60] =
            t54 * (std::pow(t104, 2) * t24 - t104 * t184 * t207 + t25 + t27);
        dA[61] = t170 * (-2 * t210 * t211 - t212);
        dA[62] = t170 * (2 * t112 * t160 * t211 - t214);
        dA[63] = t54
            * (-t109 + t110 * t24 * t42 - t17 + 2 * t209 * t24 * t42
               - t209 * t45);
        dA[64] = t54 * (-t115 * t209 + t117 * t215 + t139);
        dA[65] = t54 * (t110 * t75 - t142 * t209 - t155 + t21 + t215 * t75);
        dA[66] = t54 * (-t110 * t45 - t186 - t215 * t82);
        dA[67] = -t54 * (t110 * t115 + t197 + t215 * t89);
        dA[68] = t54 * (-t110 * t142 - t204 - t215 * t97);
        dA[69] = t170 * (-t212 - 2 * t213 * t216);
        dA[70] =
            t54 * (std::pow(t110, 2) * t24 + t110 * t184 * t209 + t25 + t26);
        dA[71] = t170 * (2 * t112 * t160 * t216 - t217);
        dA[72] = t54 * (-t111 - t113 - t12 - t218);
        dA[73] = t54 * (-t140 - t141 - t21 - t219);
        dA[74] = t54 * (-t157 + t158 + t220);
        dA[75] = t54 * (t11 + t187 + t218 + 2 * t6);
        dA[76] = t54 * (t198 + t199 + t219);
        dA[77] = t54 * (t205 - t206 - t220);
        dA[78] = -t170 * (t213 * t221 + t214);
        dA[79] = -t170 * (t210 * t221 + t217);
        dA[80] = t54 * (-std::pow(t112, 2) * t24 + t28);
    }
}
