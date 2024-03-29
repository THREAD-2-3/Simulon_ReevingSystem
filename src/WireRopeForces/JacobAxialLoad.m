function Pq = JacobAxialLoad(q,s,EA,reevSys)

rx1 = q(1);  ry1  = q(2); rz1 = q(3); rx2 = q(4); ry2 = q(5); rz2  = q(6);
qx1 = q(7); qx2  = q(8); qx3  = q(9); qy1  = q(10); qy2  = q(11); qy3  = q(12); qz1  = q(13); qz2  = q(14); qz3  = q(15);
s1 = q(16); s2 = q(17);

t2 = pi.^2;
t3 = rx1.*2.0;
t4 = rx2.*2.0;
t5 = ry1.*2.0;
t6 = ry2.*2.0;
t7 = rz1.*2.0;
t8 = rz2.*2.0;
t9 = -rx2;
t11 = -ry2;
t13 = -rz2;
t15 = -s1;
t16 = -s2;
t10 = -t4;
t12 = -t6;
t14 = -t8;
t17 = rx1+t9;
t18 = ry1+t11;
t19 = rz1+t13;
t20 = s+t15;
t21 = s1+t16;
t22 = t3+t10;
t23 = t5+t12;
t24 = t17.^2;
t25 = t7+t14;
t26 = t18.^2;
t27 = t19.^2;
t28 = 1.0./t21;
t29 = t28.^2;
t30 = t28.^3;
t31 = t28.*pi;
t32 = rx1.*t28;
t34 = rx2.*t28;
t36 = ry1.*t28;
t38 = ry2.*t28;
t40 = rz1.*t28;
t42 = rz2.*t28;
t46 = t9.*t28;
t48 = t11.*t28;
t50 = t13.*t28;
t64 = t24+t26+t27;
t33 = rx1.*t29;
t35 = rx2.*t29;
t37 = ry1.*t29;
t39 = ry2.*t29;
t41 = rz1.*t29;
t43 = rz2.*t29;
t44 = t31.*2.0;
t45 = t31.*3.0;
t47 = t9.*t29;
t49 = t11.*t29;
t51 = t13.*t29;
t52 = t20.*t31;
t53 = t20.*t29.*pi;
t66 = 1.0./sqrt(t64);
t54 = t20.*t44;
t55 = t20.*t45;
t56 = t53.*2.0;
t57 = t53.*3.0;
t58 = cos(t52);
t59 = sin(t52);
t65 = t31+t53;
t67 = t66.^3;
t60 = cos(t54);
t61 = cos(t55);
t62 = sin(t54);
t63 = sin(t55);
t68 = t44+t56;
t69 = t45+t57;
t70 = qx1.*t31.*t58.*t66;
t74 = qx1.*t17.*t29.*t58.*t66.*pi;
t76 = qx1.*t18.*t29.*t58.*t66.*pi;
t78 = qx1.*t19.*t29.*t58.*t66.*pi;
t100 = (qx1.*t17.*t22.*t31.*t58.*t67)./2.0;
t102 = (qx1.*t18.*t22.*t31.*t58.*t67)./2.0;
t103 = (qx1.*t17.*t23.*t31.*t58.*t67)./2.0;
t104 = (qx1.*t19.*t22.*t31.*t58.*t67)./2.0;
t105 = (qx1.*t18.*t23.*t31.*t58.*t67)./2.0;
t106 = (qx1.*t17.*t25.*t31.*t58.*t67)./2.0;
t108 = (qx1.*t19.*t23.*t31.*t58.*t67)./2.0;
t109 = (qx1.*t18.*t25.*t31.*t58.*t67)./2.0;
t110 = (qx1.*t19.*t25.*t31.*t58.*t67)./2.0;
t71 = qx2.*t44.*t60.*t66;
t72 = qx3.*t45.*t61.*t66;
t73 = t17.*t70;
t75 = t18.*t70;
t77 = t19.*t70;
t80 = qx2.*t17.*t29.*t60.*t66.*pi.*2.0;
t82 = qx3.*t17.*t29.*t61.*t66.*pi.*3.0;
t84 = qx2.*t18.*t29.*t60.*t66.*pi.*2.0;
t86 = qx3.*t18.*t29.*t61.*t66.*pi.*3.0;
t88 = qx2.*t19.*t29.*t60.*t66.*pi.*2.0;
t90 = qx3.*t19.*t29.*t61.*t66.*pi.*3.0;
t91 = qx2.*t17.*t22.*t31.*t60.*t67;
t92 = qx2.*t18.*t22.*t31.*t60.*t67;
t93 = qx2.*t17.*t23.*t31.*t60.*t67;
t94 = qx2.*t19.*t22.*t31.*t60.*t67;
t95 = qx2.*t18.*t23.*t31.*t60.*t67;
t96 = qx2.*t17.*t25.*t31.*t60.*t67;
t97 = qx2.*t19.*t23.*t31.*t60.*t67;
t98 = qx2.*t18.*t25.*t31.*t60.*t67;
t99 = qx2.*t19.*t25.*t31.*t60.*t67;
t112 = -t100;
t113 = qx3.*t17.*t22.*t31.*t61.*t67.*(3.0./2.0);
t114 = qx3.*t18.*t22.*t31.*t61.*t67.*(3.0./2.0);
t115 = qx3.*t17.*t23.*t31.*t61.*t67.*(3.0./2.0);
t116 = -t105;
t117 = qx3.*t19.*t22.*t31.*t61.*t67.*(3.0./2.0);
t118 = qx3.*t18.*t23.*t31.*t61.*t67.*(3.0./2.0);
t119 = qx3.*t17.*t25.*t31.*t61.*t67.*(3.0./2.0);
t120 = qx3.*t19.*t23.*t31.*t61.*t67.*(3.0./2.0);
t121 = qx3.*t18.*t25.*t31.*t61.*t67.*(3.0./2.0);
t122 = -t110;
t123 = qx3.*t19.*t25.*t31.*t61.*t67.*(3.0./2.0);
t79 = t17.*t71;
t81 = t17.*t72;
t83 = t18.*t71;
t85 = t18.*t72;
t87 = t19.*t71;
t89 = t19.*t72;
t101 = -t91;
t107 = -t95;
t111 = -t99;
t124 = -t113;
t125 = -t118;
t126 = -t123;
t130 = t92+t102+t114;
t131 = t93+t103+t115;
t132 = t94+t104+t117;
t133 = t96+t106+t119;
t134 = t97+t108+t120;
t135 = t98+t109+t121;
t127 = t32+t46+t73+t79+t81;
t128 = t36+t48+t75+t83+t85;
t129 = t40+t50+t77+t87+t89;
t136 = t28+t70+t71+t72+t101+t112+t124;
t137 = t28+t70+t71+t72+t107+t116+t125;
t138 = t28+t70+t71+t72+t111+t122+t126;
t139 = t127.*t131;
t140 = t127.*t133;
t141 = t128.*t130;
t142 = t128.*t135;
t143 = t129.*t132;
t144 = t129.*t134;
t145 = t127.*t136;
t146 = t128.*t137;
t147 = t129.*t138;
t148 = -t145;
t149 = -t146;
t150 = -t147;
t151 = t141+t143+t148;
t152 = t139+t144+t149;
t153 = t140+t142+t150;
t154 = EA.*t151;
t155 = EA.*t152;
t156 = EA.*t153;
Pq = [-t154,-t155,-t156,t154,t155,t156,EA.*(t17.*t31.*t58.*t66.*t127+t18.*t31.*t58.*t66.*t128+t19.*t31.*t58.*t66.*t129),EA.*(t17.*t44.*t60.*t66.*t127+t18.*t44.*t60.*t66.*t128+t19.*t44.*t60.*t66.*t129),EA.*(t17.*t45.*t61.*t66.*t127+t18.*t45.*t61.*t66.*t128+t19.*t45.*t61.*t66.*t129),0.0,0.0,0.0,0.0,0.0,0.0,-EA.*(t127.*(t33+t47+t74+t80+t82-qx1.*t17.*t31.*t59.*t65.*t66-qx2.*t17.*t31.*t62.*t66.*t68.*2.0-qx3.*t17.*t31.*t63.*t66.*t69.*3.0)+t128.*(t37+t49+t76+t84+t86-qx1.*t18.*t31.*t59.*t65.*t66-qx2.*t18.*t31.*t62.*t66.*t68.*2.0-qx3.*t18.*t31.*t63.*t66.*t69.*3.0)+t129.*(t41+t51+t78+t88+t90-qx1.*t19.*t31.*t59.*t65.*t66-qx2.*t19.*t31.*t62.*t66.*t68.*2.0-qx3.*t19.*t31.*t63.*t66.*t69.*3.0)),EA.*(t127.*(t33+t47+t74+t80+t82-qx1.*t2.*t17.*t20.*t30.*t59.*t66-qx2.*t2.*t17.*t20.*t30.*t62.*t66.*4.0-qx3.*t2.*t17.*t20.*t30.*t63.*t66.*9.0)+t128.*(t37+t49+t76+t84+t86-qx1.*t2.*t18.*t20.*t30.*t59.*t66-qx2.*t2.*t18.*t20.*t30.*t62.*t66.*4.0-qx3.*t2.*t18.*t20.*t30.*t63.*t66.*9.0)+t129.*(t41+t51+t78+t88+t90-qx1.*t2.*t19.*t20.*t30.*t59.*t66-qx2.*t2.*t19.*t20.*t30.*t62.*t66.*4.0-qx3.*t2.*t19.*t20.*t30.*t63.*t66.*9.0))];
