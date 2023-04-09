function Dudef = jacobDiffDefEnergy_noDamp(q,s,EA,EI)


qq = [q' s]';

Dudef = MexjacobDiffDefEnergy(qq,EI,EA);


% 
% rx1 = q(1);  ry1  = q(2); rz1 = q(3); rx2 = q(4); ry2 = q(5); rz2  = q(6);
% qx1 = q(7); qx2  = q(8); qx3  = q(9); qy1  = q(10); qy2  = q(11); qy3  = q(12); qz1  = q(13); qz2  = q(14); qz3  = q(15);
% s1 = q(16); s2 = q(17);
% 
% t2 = pi.^2;
% t3 = pi.^3;
% t5 = pi.^5;
% t6 = rx1.*2.0;
% t7 = rx2.*2.0;
% t8 = ry1.*2.0;
% t9 = ry2.*2.0;
% t10 = rz1.*2.0;
% t11 = rz2.*2.0;
% t12 = -rx2;
% t14 = -ry2;
% t16 = -rz2;
% t18 = -s1;
% t19 = -s2;
% t4 = t2.^2;
% t13 = -t7;
% t15 = -t9;
% t17 = -t11;
% t20 = rx1+t12;
% t21 = ry1+t14;
% t22 = rz1+t16;
% t23 = s+t18;
% t24 = s1+t19;
% t25 = t6+t13;
% t26 = t8+t15;
% t27 = t20.^2;
% t28 = t10+t17;
% t29 = t21.^2;
% t30 = t22.^2;
% t31 = 1.0./t24;
% t32 = t31.^2;
% t33 = t31.^3;
% t35 = t31.^5;
% t37 = t31.*pi;
% t38 = rx1.*t31;
% t40 = rx2.*t31;
% t42 = ry1.*t31;
% t44 = ry2.*t31;
% t46 = rz1.*t31;
% t48 = rz2.*t31;
% t52 = t12.*t31;
% t54 = t14.*t31;
% t56 = t16.*t31;
% t76 = t27+t29+t30;
% t34 = t32.^2;
% t36 = t32.^3;
% t39 = rx1.*t32;
% t41 = rx2.*t32;
% t43 = ry1.*t32;
% t45 = ry2.*t32;
% t47 = rz1.*t32;
% t49 = rz2.*t32;
% t50 = t37.*2.0;
% t51 = t37.*3.0;
% t53 = t12.*t32;
% t55 = t14.*t32;
% t57 = t16.*t32;
% t58 = t23.*t37;
% t59 = t23.*t32.*pi;
% t78 = 1.0./sqrt(t76);
% t60 = t23.*t50;
% t61 = t23.*t51;
% t62 = t59.*2.0;
% t63 = t59.*3.0;
% t64 = cos(t58);
% t65 = sin(t58);
% t77 = t37+t59;
% t79 = t78.^3;
% t66 = cos(t60);
% t67 = cos(t61);
% t68 = sin(t60);
% t69 = sin(t61);
% t70 = t64.^2;
% t71 = t65.^2;
% t80 = t50+t62;
% t81 = t51+t63;
% t100 = qx1.*t37.*t64.*t78;
% t104 = qx1.*t20.*t32.*t64.*t78.*pi;
% t106 = qx1.*t21.*t32.*t64.*t78.*pi;
% t108 = qx1.*t22.*t32.*t64.*t78.*pi;
% t121 = qx1.*t2.*t20.*t23.*t33.*t65.*t78;
% t122 = qx1.*t2.*t21.*t23.*t33.*t65.*t78;
% t123 = qx1.*t2.*t22.*t23.*t33.*t65.*t78;
% t139 = (qx1.*t20.*t25.*t37.*t64.*t79)./2.0;
% t141 = (qx1.*t21.*t25.*t37.*t64.*t79)./2.0;
% t142 = (qx1.*t20.*t26.*t37.*t64.*t79)./2.0;
% t143 = (qx1.*t22.*t25.*t37.*t64.*t79)./2.0;
% t144 = (qx1.*t21.*t26.*t37.*t64.*t79)./2.0;
% t145 = (qx1.*t20.*t28.*t37.*t64.*t79)./2.0;
% t147 = (qx1.*t22.*t26.*t37.*t64.*t79)./2.0;
% t148 = (qx1.*t21.*t28.*t37.*t64.*t79)./2.0;
% t149 = (qx1.*t22.*t28.*t37.*t64.*t79)./2.0;
% t175 = qx1.*t20.*t37.*t65.*t77.*t78;
% t176 = qx1.*t21.*t37.*t65.*t77.*t78;
% t177 = qx1.*t22.*t37.*t65.*t77.*t78;
% t72 = t66.^2;
% t73 = t67.^2;
% t74 = t68.^2;
% t75 = t69.^2;
% t82 = EI.*qy1.*t4.*t35.*t71.*2.0;
% t83 = EI.*qz1.*t4.*t35.*t71.*2.0;
% t88 = EI.*qy1.*t4.*t35.*t65.*t68.*8.0;
% t89 = EI.*qy2.*t4.*t35.*t65.*t68.*8.0;
% t90 = EI.*qy1.*t4.*t35.*t65.*t69.*1.8e+1;
% t91 = EI.*qy3.*t4.*t35.*t65.*t69.*1.8e+1;
% t92 = EI.*qz1.*t4.*t35.*t65.*t68.*8.0;
% t93 = EI.*qz2.*t4.*t35.*t65.*t68.*8.0;
% t94 = EI.*qz1.*t4.*t35.*t65.*t69.*1.8e+1;
% t95 = EI.*qz3.*t4.*t35.*t65.*t69.*1.8e+1;
% t96 = EI.*qy2.*t4.*t35.*t68.*t69.*7.2e+1;
% t97 = EI.*qy3.*t4.*t35.*t68.*t69.*7.2e+1;
% t98 = EI.*qz2.*t4.*t35.*t68.*t69.*7.2e+1;
% t99 = EI.*qz3.*t4.*t35.*t68.*t69.*7.2e+1;
% t101 = qx2.*t50.*t66.*t78;
% t102 = qx3.*t51.*t67.*t78;
% t103 = t20.*t100;
% t105 = t21.*t100;
% t107 = t22.*t100;
% t110 = qx2.*t20.*t32.*t66.*t78.*pi.*2.0;
% t112 = qx3.*t20.*t32.*t67.*t78.*pi.*3.0;
% t114 = qx2.*t21.*t32.*t66.*t78.*pi.*2.0;
% t116 = qx3.*t21.*t32.*t67.*t78.*pi.*3.0;
% t118 = qx2.*t22.*t32.*t66.*t78.*pi.*2.0;
% t120 = qx3.*t22.*t32.*t67.*t78.*pi.*3.0;
% t124 = qx2.*t20.*t25.*t37.*t66.*t79;
% t125 = qx2.*t21.*t25.*t37.*t66.*t79;
% t126 = qx2.*t20.*t26.*t37.*t66.*t79;
% t127 = qx2.*t22.*t25.*t37.*t66.*t79;
% t128 = qx2.*t21.*t26.*t37.*t66.*t79;
% t129 = qx2.*t20.*t28.*t37.*t66.*t79;
% t130 = qx2.*t22.*t26.*t37.*t66.*t79;
% t131 = qx2.*t21.*t28.*t37.*t66.*t79;
% t132 = qx2.*t22.*t28.*t37.*t66.*t79;
% t133 = -t121;
% t134 = qx2.*t2.*t20.*t23.*t33.*t68.*t78.*4.0;
% t135 = -t122;
% t136 = qx2.*t2.*t21.*t23.*t33.*t68.*t78.*4.0;
% t137 = -t123;
% t138 = qx2.*t2.*t22.*t23.*t33.*t68.*t78.*4.0;
% t152 = qx3.*t2.*t20.*t23.*t33.*t69.*t78.*9.0;
% t154 = qx3.*t2.*t21.*t23.*t33.*t69.*t78.*9.0;
% t156 = qx3.*t2.*t22.*t23.*t33.*t69.*t78.*9.0;
% t157 = -t139;
% t158 = qx3.*t20.*t25.*t37.*t67.*t79.*(3.0./2.0);
% t159 = qx3.*t21.*t25.*t37.*t67.*t79.*(3.0./2.0);
% t160 = qx3.*t20.*t26.*t37.*t67.*t79.*(3.0./2.0);
% t161 = -t144;
% t162 = qx3.*t22.*t25.*t37.*t67.*t79.*(3.0./2.0);
% t163 = qx3.*t21.*t26.*t37.*t67.*t79.*(3.0./2.0);
% t164 = qx3.*t20.*t28.*t37.*t67.*t79.*(3.0./2.0);
% t165 = qx3.*t22.*t26.*t37.*t67.*t79.*(3.0./2.0);
% t166 = qx3.*t21.*t28.*t37.*t67.*t79.*(3.0./2.0);
% t167 = -t149;
% t168 = qx3.*t22.*t28.*t37.*t67.*t79.*(3.0./2.0);
% t178 = -t175;
% t179 = -t176;
% t180 = -t177;
% t181 = qx2.*t20.*t50.*t68.*t78.*t80;
% t182 = qx3.*t20.*t51.*t69.*t78.*t81;
% t183 = qx2.*t21.*t50.*t68.*t78.*t80;
% t184 = qx3.*t21.*t51.*t69.*t78.*t81;
% t185 = qx2.*t22.*t50.*t68.*t78.*t80;
% t186 = qx3.*t22.*t51.*t69.*t78.*t81;
% t187 = qx2.*t20.*t37.*t68.*t78.*t80.*-2.0;
% t188 = qx3.*t20.*t37.*t69.*t78.*t81.*-3.0;
% t189 = qx2.*t21.*t37.*t68.*t78.*t80.*-2.0;
% t190 = qx3.*t21.*t37.*t69.*t78.*t81.*-3.0;
% t191 = qx2.*t22.*t37.*t68.*t78.*t80.*-2.0;
% t192 = qx3.*t22.*t37.*t69.*t78.*t81.*-3.0;
% t84 = EI.*qy2.*t4.*t35.*t74.*3.2e+1;
% t85 = EI.*qz2.*t4.*t35.*t74.*3.2e+1;
% t86 = EI.*qy3.*t4.*t35.*t75.*1.62e+2;
% t87 = EI.*qz3.*t4.*t35.*t75.*1.62e+2;
% t109 = t20.*t101;
% t111 = t20.*t102;
% t113 = t21.*t101;
% t115 = t21.*t102;
% t117 = t22.*t101;
% t119 = t22.*t102;
% t140 = -t124;
% t146 = -t128;
% t150 = -t132;
% t151 = -t134;
% t153 = -t136;
% t155 = -t138;
% t169 = -t152;
% t170 = -t154;
% t171 = -t156;
% t172 = -t158;
% t173 = -t163;
% t174 = -t168;
% t199 = t125+t141+t159;
% t200 = t126+t142+t160;
% t201 = t127+t143+t162;
% t202 = t129+t145+t164;
% t203 = t130+t147+t165;
% t204 = t131+t148+t166;
% t229 = t39+t53+t104+t110+t112+t178+t187+t188;
% t230 = t43+t55+t106+t114+t116+t179+t189+t190;
% t231 = t47+t57+t108+t118+t120+t180+t191+t192;
% t193 = t38+t52+t103+t109+t111;
% t194 = t42+t54+t105+t113+t115;
% t195 = t46+t56+t107+t117+t119;
% t217 = t31+t100+t101+t102+t140+t157+t172;
% t218 = t31+t100+t101+t102+t146+t161+t173;
% t219 = t31+t100+t101+t102+t150+t167+t174;
% t220 = t39+t53+t104+t110+t112+t133+t151+t169;
% t221 = t43+t55+t106+t114+t116+t135+t153+t170;
% t222 = t47+t57+t108+t118+t120+t137+t155+t171;
% t196 = t193.^2;
% t197 = t194.^2;
% t198 = t195.^2;
% t208 = t20.*t37.*t64.*t78.*t193;
% t209 = t21.*t37.*t64.*t78.*t194;
% t210 = t22.*t37.*t64.*t78.*t195;
% t211 = t20.*t50.*t66.*t78.*t193;
% t212 = t20.*t51.*t67.*t78.*t193;
% t213 = t21.*t50.*t66.*t78.*t194;
% t214 = t21.*t51.*t67.*t78.*t194;
% t215 = t22.*t50.*t66.*t78.*t195;
% t216 = t22.*t51.*t67.*t78.*t195;
% t223 = t193.*t200;
% t224 = t193.*t202;
% t225 = t194.*t199;
% t226 = t194.*t204;
% t227 = t195.*t201;
% t228 = t195.*t203;
% t232 = t193.*t217;
% t233 = t194.*t218;
% t234 = t195.*t219;
% t238 = t193.*t220;
% t239 = t194.*t221;
% t240 = t195.*t222;
% t248 = t193.*t229;
% t249 = t194.*t230;
% t250 = t195.*t231;
% t205 = t196./2.0;
% t206 = t197./2.0;
% t207 = t198./2.0;
% t235 = -t232;
% t236 = -t233;
% t237 = -t234;
% t263 = t208+t209+t210;
% t264 = t211+t213+t215;
% t265 = t212+t214+t216;
% t323 = t238+t239+t240;
% t324 = t248+t249+t250;
% t241 = t205+t206+t207-1.0./2.0;
% t266 = t225+t227+t235;
% t267 = t223+t228+t236;
% t268 = t224+t226+t237;
% t242 = EA.*qy1.*t2.*t33.*t70.*t241;
% t243 = EA.*qz1.*t2.*t33.*t70.*t241;
% t244 = EA.*qy2.*t2.*t33.*t72.*t241.*4.0;
% t245 = EA.*qz2.*t2.*t33.*t72.*t241.*4.0;
% t246 = EA.*qy3.*t2.*t33.*t73.*t241.*9.0;
% t247 = EA.*qz3.*t2.*t33.*t73.*t241.*9.0;
% t251 = EA.*qy1.*t2.*t33.*t64.*t66.*t241.*2.0;
% t252 = EA.*qy2.*t2.*t33.*t64.*t66.*t241.*2.0;
% t253 = EA.*qy1.*t2.*t33.*t64.*t67.*t241.*3.0;
% t254 = EA.*qy3.*t2.*t33.*t64.*t67.*t241.*3.0;
% t255 = EA.*qz1.*t2.*t33.*t64.*t66.*t241.*2.0;
% t256 = EA.*qz2.*t2.*t33.*t64.*t66.*t241.*2.0;
% t257 = EA.*qz1.*t2.*t33.*t64.*t67.*t241.*3.0;
% t258 = EA.*qz3.*t2.*t33.*t64.*t67.*t241.*3.0;
% t259 = EA.*qy2.*t2.*t33.*t66.*t67.*t241.*6.0;
% t260 = EA.*qy3.*t2.*t33.*t66.*t67.*t241.*6.0;
% t261 = EA.*qz2.*t2.*t33.*t66.*t67.*t241.*6.0;
% t262 = EA.*qz3.*t2.*t33.*t66.*t67.*t241.*6.0;
% t269 = EA.*qy2.*t2.*t32.*t72.*t266.*2.0;
% t270 = EA.*qz2.*t2.*t32.*t72.*t266.*2.0;
% t271 = EA.*qy2.*t2.*t32.*t72.*t267.*2.0;
% t272 = EA.*qz2.*t2.*t32.*t72.*t267.*2.0;
% t273 = EA.*qy2.*t2.*t32.*t72.*t268.*2.0;
% t274 = EA.*qz2.*t2.*t32.*t72.*t268.*2.0;
% t275 = (EA.*qy1.*t2.*t32.*t70.*t266)./2.0;
% t276 = (EA.*qz1.*t2.*t32.*t70.*t266)./2.0;
% t277 = (EA.*qy1.*t2.*t32.*t70.*t267)./2.0;
% t278 = (EA.*qz1.*t2.*t32.*t70.*t267)./2.0;
% t279 = (EA.*qy1.*t2.*t32.*t70.*t268)./2.0;
% t280 = (EA.*qz1.*t2.*t32.*t70.*t268)./2.0;
% t281 = EA.*qy3.*t2.*t32.*t73.*t266.*(9.0./2.0);
% t282 = EA.*qz3.*t2.*t32.*t73.*t266.*(9.0./2.0);
% t283 = EA.*qy3.*t2.*t32.*t73.*t267.*(9.0./2.0);
% t284 = EA.*qz3.*t2.*t32.*t73.*t267.*(9.0./2.0);
% t285 = EA.*qy3.*t2.*t32.*t73.*t268.*(9.0./2.0);
% t286 = EA.*qz3.*t2.*t32.*t73.*t268.*(9.0./2.0);
% t287 = EA.*qy1.*t2.*t32.*t64.*t66.*t266;
% t288 = EA.*qy2.*t2.*t32.*t64.*t66.*t266;
% t289 = EA.*qz1.*t2.*t32.*t64.*t66.*t266;
% t290 = EA.*qz2.*t2.*t32.*t64.*t66.*t266;
% t291 = EA.*qy1.*t2.*t32.*t64.*t66.*t267;
% t292 = EA.*qy2.*t2.*t32.*t64.*t66.*t267;
% t293 = EA.*qz1.*t2.*t32.*t64.*t66.*t267;
% t294 = EA.*qz2.*t2.*t32.*t64.*t66.*t267;
% t295 = EA.*qy1.*t2.*t32.*t64.*t66.*t268;
% t296 = EA.*qy2.*t2.*t32.*t64.*t66.*t268;
% t297 = EA.*qz1.*t2.*t32.*t64.*t66.*t268;
% t298 = EA.*qz2.*t2.*t32.*t64.*t66.*t268;
% t299 = EA.*qy2.*t2.*t32.*t66.*t67.*t266.*3.0;
% t300 = EA.*qy3.*t2.*t32.*t66.*t67.*t266.*3.0;
% t301 = EA.*qz2.*t2.*t32.*t66.*t67.*t266.*3.0;
% t302 = EA.*qz3.*t2.*t32.*t66.*t67.*t266.*3.0;
% t303 = EA.*qy2.*t2.*t32.*t66.*t67.*t267.*3.0;
% t304 = EA.*qy3.*t2.*t32.*t66.*t67.*t267.*3.0;
% t305 = EA.*qz2.*t2.*t32.*t66.*t67.*t267.*3.0;
% t306 = EA.*qz3.*t2.*t32.*t66.*t67.*t267.*3.0;
% t307 = EA.*qy2.*t2.*t32.*t66.*t67.*t268.*3.0;
% t308 = EA.*qy3.*t2.*t32.*t66.*t67.*t268.*3.0;
% t309 = EA.*qz2.*t2.*t32.*t66.*t67.*t268.*3.0;
% t310 = EA.*qz3.*t2.*t32.*t66.*t67.*t268.*3.0;
% t311 = EA.*qy1.*t2.*t32.*t64.*t67.*t266.*(3.0./2.0);
% t312 = EA.*qy3.*t2.*t32.*t64.*t67.*t266.*(3.0./2.0);
% t313 = EA.*qz1.*t2.*t32.*t64.*t67.*t266.*(3.0./2.0);
% t314 = EA.*qz3.*t2.*t32.*t64.*t67.*t266.*(3.0./2.0);
% t315 = EA.*qy1.*t2.*t32.*t64.*t67.*t267.*(3.0./2.0);
% t316 = EA.*qy3.*t2.*t32.*t64.*t67.*t267.*(3.0./2.0);
% t317 = EA.*qz1.*t2.*t32.*t64.*t67.*t267.*(3.0./2.0);
% t318 = EA.*qz3.*t2.*t32.*t64.*t67.*t267.*(3.0./2.0);
% t319 = EA.*qy1.*t2.*t32.*t64.*t67.*t268.*(3.0./2.0);
% t320 = EA.*qy3.*t2.*t32.*t64.*t67.*t268.*(3.0./2.0);
% t321 = EA.*qz1.*t2.*t32.*t64.*t67.*t268.*(3.0./2.0);
% t322 = EA.*qz3.*t2.*t32.*t64.*t67.*t268.*(3.0./2.0);
% t325 = EA.*t241.*t266;
% t326 = EA.*t241.*t267;
% t327 = EA.*t241.*t268;
% t328 = t269+t287+t300;
% t329 = t270+t289+t302;
% t330 = t271+t291+t304;
% t331 = t272+t293+t306;
% t332 = t273+t295+t308;
% t333 = t274+t297+t310;
% t340 = t275+t288+t312;
% t341 = t276+t290+t314;
% t342 = t277+t292+t316;
% t343 = t278+t294+t318;
% t344 = t279+t296+t320;
% t345 = t280+t298+t322;
% t352 = t281+t299+t311;
% t353 = t282+t301+t313;
% t354 = t283+t303+t315;
% t355 = t284+t305+t317;
% t356 = t285+t307+t319;
% t357 = t286+t309+t321;
% t334 = qy2.*t328;
% t335 = qz2.*t329;
% t336 = qy2.*t330;
% t337 = qz2.*t331;
% t338 = qy2.*t332;
% t339 = qz2.*t333;
% t346 = qy1.*t340;
% t347 = qz1.*t341;
% t348 = qy1.*t342;
% t349 = qz1.*t343;
% t350 = qy1.*t344;
% t351 = qz1.*t345;
% t358 = qy3.*t352;
% t359 = qz3.*t353;
% t360 = qy3.*t354;
% t361 = qz3.*t355;
% t362 = qy3.*t356;
% t363 = qz3.*t357;
% Dudef = [-t325-t334-t335-t346-t347-t358-t359,-t326-t336-t337-t348-t349-t360-t361,-t327-t338-t339-t350-t351-t362-t363,t325+t334+t335+t346+t347+t358+t359,t326+t336+t337+t348+t349+t360+t361,t327+t338+t339+t350+t351+t362+t363,qy1.*((EA.*qy1.*t2.*t32.*t70.*t263)./2.0+EA.*qy2.*t2.*t32.*t64.*t66.*t263+EA.*qy3.*t2.*t32.*t64.*t67.*t263.*(3.0./2.0))+qy2.*(EA.*qy2.*t2.*t32.*t72.*t263.*2.0+EA.*qy1.*t2.*t32.*t64.*t66.*t263+EA.*qy3.*t2.*t32.*t66.*t67.*t263.*3.0)+qy3.*(EA.*qy3.*t2.*t32.*t73.*t263.*(9.0./2.0)+EA.*qy1.*t2.*t32.*t64.*t67.*t263.*(3.0./2.0)+EA.*qy2.*t2.*t32.*t66.*t67.*t263.*3.0)+qz1.*((EA.*qz1.*t2.*t32.*t70.*t263)./2.0+EA.*qz2.*t2.*t32.*t64.*t66.*t263+EA.*qz3.*t2.*t32.*t64.*t67.*t263.*(3.0./2.0))+qz2.*(EA.*qz2.*t2.*t32.*t72.*t263.*2.0+EA.*qz1.*t2.*t32.*t64.*t66.*t263+EA.*qz3.*t2.*t32.*t66.*t67.*t263.*3.0)+qz3.*(EA.*qz3.*t2.*t32.*t73.*t263.*(9.0./2.0)+EA.*qz1.*t2.*t32.*t64.*t67.*t263.*(3.0./2.0)+EA.*qz2.*t2.*t32.*t66.*t67.*t263.*3.0)+EA.*t241.*t263,qy1.*((EA.*qy1.*t2.*t32.*t70.*t264)./2.0+EA.*qy2.*t2.*t32.*t64.*t66.*t264+EA.*qy3.*t2.*t32.*t64.*t67.*t264.*(3.0./2.0))+qy2.*(EA.*qy2.*t2.*t32.*t72.*t264.*2.0+EA.*qy1.*t2.*t32.*t64.*t66.*t264+EA.*qy3.*t2.*t32.*t66.*t67.*t264.*3.0)+qy3.*(EA.*qy3.*t2.*t32.*t73.*t264.*(9.0./2.0)+EA.*qy1.*t2.*t32.*t64.*t67.*t264.*(3.0./2.0)+EA.*qy2.*t2.*t32.*t66.*t67.*t264.*3.0)+qz1.*((EA.*qz1.*t2.*t32.*t70.*t264)./2.0+EA.*qz2.*t2.*t32.*t64.*t66.*t264+EA.*qz3.*t2.*t32.*t64.*t67.*t264.*(3.0./2.0))+qz2.*(EA.*qz2.*t2.*t32.*t72.*t264.*2.0+EA.*qz1.*t2.*t32.*t64.*t66.*t264+EA.*qz3.*t2.*t32.*t66.*t67.*t264.*3.0)+qz3.*(EA.*qz3.*t2.*t32.*t73.*t264.*(9.0./2.0)+EA.*qz1.*t2.*t32.*t64.*t67.*t264.*(3.0./2.0)+EA.*qz2.*t2.*t32.*t66.*t67.*t264.*3.0)+EA.*t241.*t264,qy1.*((EA.*qy1.*t2.*t32.*t70.*t265)./2.0+EA.*qy2.*t2.*t32.*t64.*t66.*t265+EA.*qy3.*t2.*t32.*t64.*t67.*t265.*(3.0./2.0))+qy2.*(EA.*qy2.*t2.*t32.*t72.*t265.*2.0+EA.*qy1.*t2.*t32.*t64.*t66.*t265+EA.*qy3.*t2.*t32.*t66.*t67.*t265.*3.0)+qy3.*(EA.*qy3.*t2.*t32.*t73.*t265.*(9.0./2.0)+EA.*qy1.*t2.*t32.*t64.*t67.*t265.*(3.0./2.0)+EA.*qy2.*t2.*t32.*t66.*t67.*t265.*3.0)+qz1.*((EA.*qz1.*t2.*t32.*t70.*t265)./2.0+EA.*qz2.*t2.*t32.*t64.*t66.*t265+EA.*qz3.*t2.*t32.*t64.*t67.*t265.*(3.0./2.0))+qz2.*(EA.*qz2.*t2.*t32.*t72.*t265.*2.0+EA.*qz1.*t2.*t32.*t64.*t66.*t265+EA.*qz3.*t2.*t32.*t66.*t67.*t265.*3.0)+qz3.*(EA.*qz3.*t2.*t32.*t73.*t265.*(9.0./2.0)+EA.*qz1.*t2.*t32.*t64.*t67.*t265.*(3.0./2.0)+EA.*qz2.*t2.*t32.*t66.*t67.*t265.*3.0)+EA.*t241.*t265,EI.*qy1.*t4.*t34.*t71+EA.*qy1.*t2.*t32.*t70.*t241+EI.*qy2.*t4.*t34.*t65.*t68.*4.0+EI.*qy3.*t4.*t34.*t65.*t69.*9.0+EA.*qy2.*t2.*t32.*t64.*t66.*t241.*2.0+EA.*qy3.*t2.*t32.*t64.*t67.*t241.*3.0,EI.*qy2.*t4.*t34.*t74.*1.6e+1+EA.*qy2.*t2.*t32.*t72.*t241.*4.0+EI.*qy1.*t4.*t34.*t65.*t68.*4.0+EI.*qy3.*t4.*t34.*t68.*t69.*3.6e+1+EA.*qy1.*t2.*t32.*t64.*t66.*t241.*2.0+EA.*qy3.*t2.*t32.*t66.*t67.*t241.*6.0,EI.*qy3.*t4.*t34.*t75.*8.1e+1+EA.*qy3.*t2.*t32.*t73.*t241.*9.0+EI.*qy1.*t4.*t34.*t65.*t69.*9.0+EI.*qy2.*t4.*t34.*t68.*t69.*3.6e+1+EA.*qy1.*t2.*t32.*t64.*t67.*t241.*3.0+EA.*qy2.*t2.*t32.*t66.*t67.*t241.*6.0,EI.*qz1.*t4.*t34.*t71+EA.*qz1.*t2.*t32.*t70.*t241+EI.*qz2.*t4.*t34.*t65.*t68.*4.0+EI.*qz3.*t4.*t34.*t65.*t69.*9.0+EA.*qz2.*t2.*t32.*t64.*t66.*t241.*2.0+EA.*qz3.*t2.*t32.*t64.*t67.*t241.*3.0,EI.*qz2.*t4.*t34.*t74.*1.6e+1+EA.*qz2.*t2.*t32.*t72.*t241.*4.0+EI.*qz1.*t4.*t34.*t65.*t68.*4.0+EI.*qz3.*t4.*t34.*t68.*t69.*3.6e+1+EA.*qz1.*t2.*t32.*t64.*t66.*t241.*2.0+EA.*qz3.*t2.*t32.*t66.*t67.*t241.*6.0,EI.*qz3.*t4.*t34.*t75.*8.1e+1+EA.*qz3.*t2.*t32.*t73.*t241.*9.0+EI.*qz1.*t4.*t34.*t65.*t69.*9.0+EI.*qz2.*t4.*t34.*t68.*t69.*3.6e+1+EA.*qz1.*t2.*t32.*t64.*t67.*t241.*3.0+EA.*qz2.*t2.*t32.*t66.*t67.*t241.*6.0,-qy1.*(t242+t252+t254+(EA.*qy1.*t2.*t32.*t70.*t324)./2.0+EA.*qy2.*t2.*t32.*t64.*t66.*t324+EA.*qy3.*t2.*t32.*t64.*t67.*t324.*(3.0./2.0)-EA.*qy1.*t2.*t32.*t64.*t65.*t77.*t241-EA.*qy2.*t2.*t32.*t65.*t66.*t77.*t241-EA.*qy3.*t2.*t32.*t65.*t67.*t77.*t241.*(3.0./2.0)-EA.*qy2.*t2.*t32.*t64.*t68.*t80.*t241-EA.*qy3.*t2.*t32.*t64.*t69.*t81.*t241.*(3.0./2.0))-qy2.*(t244+t251+t260+EA.*qy2.*t2.*t32.*t72.*t324.*2.0+EA.*qy1.*t2.*t32.*t64.*t66.*t324+EA.*qy3.*t2.*t32.*t66.*t67.*t324.*3.0-EA.*qy1.*t2.*t32.*t65.*t66.*t77.*t241-EA.*qy1.*t2.*t32.*t64.*t68.*t80.*t241-EA.*qy2.*t2.*t32.*t66.*t68.*t80.*t241.*4.0-EA.*qy3.*t2.*t32.*t67.*t68.*t80.*t241.*3.0-EA.*qy3.*t2.*t32.*t66.*t69.*t81.*t241.*3.0)-qy3.*(t246+t253+t259+EA.*qy3.*t2.*t32.*t73.*t324.*(9.0./2.0)+EA.*qy1.*t2.*t32.*t64.*t67.*t324.*(3.0./2.0)+EA.*qy2.*t2.*t32.*t66.*t67.*t324.*3.0-EA.*qy1.*t2.*t32.*t65.*t67.*t77.*t241.*(3.0./2.0)-EA.*qy1.*t2.*t32.*t64.*t69.*t81.*t241.*(3.0./2.0)-EA.*qy2.*t2.*t32.*t67.*t68.*t80.*t241.*3.0-EA.*qy2.*t2.*t32.*t66.*t69.*t81.*t241.*3.0-EA.*qy3.*t2.*t32.*t67.*t69.*t81.*t241.*9.0)-qz1.*(t243+t256+t258+(EA.*qz1.*t2.*t32.*t70.*t324)./2.0+EA.*qz2.*t2.*t32.*t64.*t66.*t324+EA.*qz3.*t2.*t32.*t64.*t67.*t324.*(3.0./2.0)-EA.*qz1.*t2.*t32.*t64.*t65.*t77.*t241-EA.*qz2.*t2.*t32.*t65.*t66.*t77.*t241-EA.*qz3.*t2.*t32.*t65.*t67.*t77.*t241.*(3.0./2.0)-EA.*qz2.*t2.*t32.*t64.*t68.*t80.*t241-EA.*qz3.*t2.*t32.*t64.*t69.*t81.*t241.*(3.0./2.0))-qz2.*(t245+t255+t262+EA.*qz2.*t2.*t32.*t72.*t324.*2.0+EA.*qz1.*t2.*t32.*t64.*t66.*t324+EA.*qz3.*t2.*t32.*t66.*t67.*t324.*3.0-EA.*qz1.*t2.*t32.*t65.*t66.*t77.*t241-EA.*qz1.*t2.*t32.*t64.*t68.*t80.*t241-EA.*qz2.*t2.*t32.*t66.*t68.*t80.*t241.*4.0-EA.*qz3.*t2.*t32.*t67.*t68.*t80.*t241.*3.0-EA.*qz3.*t2.*t32.*t66.*t69.*t81.*t241.*3.0)-qz3.*(t247+t257+t261+EA.*qz3.*t2.*t32.*t73.*t324.*(9.0./2.0)+EA.*qz1.*t2.*t32.*t64.*t67.*t324.*(3.0./2.0)+EA.*qz2.*t2.*t32.*t66.*t67.*t324.*3.0-EA.*qz1.*t2.*t32.*t65.*t67.*t77.*t241.*(3.0./2.0)-EA.*qz1.*t2.*t32.*t64.*t69.*t81.*t241.*(3.0./2.0)-EA.*qz2.*t2.*t32.*t67.*t68.*t80.*t241.*3.0-EA.*qz2.*t2.*t32.*t66.*t69.*t81.*t241.*3.0-EA.*qz3.*t2.*t32.*t67.*t69.*t81.*t241.*9.0)-qy1.*(t82+t89+t91+EI.*qy1.*t4.*t34.*t64.*t65.*t77+EI.*qy2.*t4.*t34.*t64.*t68.*t77.*2.0+EI.*qy2.*t4.*t34.*t65.*t66.*t80.*2.0+EI.*qy3.*t4.*t34.*t64.*t69.*t77.*(9.0./2.0)+EI.*qy3.*t4.*t34.*t65.*t67.*t81.*(9.0./2.0))-qy2.*(t84+t88+t97+EI.*qy1.*t4.*t34.*t64.*t68.*t77.*2.0+EI.*qy1.*t4.*t34.*t65.*t66.*t80.*2.0+EI.*qy2.*t4.*t34.*t66.*t68.*t80.*1.6e+1+EI.*qy3.*t4.*t34.*t66.*t69.*t80.*1.8e+1+EI.*qy3.*t4.*t34.*t67.*t68.*t81.*1.8e+1)-qy3.*(t86+t90+t96+EI.*qy1.*t4.*t34.*t64.*t69.*t77.*(9.0./2.0)+EI.*qy1.*t4.*t34.*t65.*t67.*t81.*(9.0./2.0)+EI.*qy2.*t4.*t34.*t66.*t69.*t80.*1.8e+1+EI.*qy2.*t4.*t34.*t67.*t68.*t81.*1.8e+1+EI.*qy3.*t4.*t34.*t67.*t69.*t81.*8.1e+1)-qz1.*(t83+t93+t95+EI.*qz1.*t4.*t34.*t64.*t65.*t77+EI.*qz2.*t4.*t34.*t64.*t68.*t77.*2.0+EI.*qz2.*t4.*t34.*t65.*t66.*t80.*2.0+EI.*qz3.*t4.*t34.*t64.*t69.*t77.*(9.0./2.0)+EI.*qz3.*t4.*t34.*t65.*t67.*t81.*(9.0./2.0))-qz2.*(t85+t92+t99+EI.*qz1.*t4.*t34.*t64.*t68.*t77.*2.0+EI.*qz1.*t4.*t34.*t65.*t66.*t80.*2.0+EI.*qz2.*t4.*t34.*t66.*t68.*t80.*1.6e+1+EI.*qz3.*t4.*t34.*t66.*t69.*t80.*1.8e+1+EI.*qz3.*t4.*t34.*t67.*t68.*t81.*1.8e+1)-qz3.*(t87+t94+t98+EI.*qz1.*t4.*t34.*t64.*t69.*t77.*(9.0./2.0)+EI.*qz1.*t4.*t34.*t65.*t67.*t81.*(9.0./2.0)+EI.*qz2.*t4.*t34.*t66.*t69.*t80.*1.8e+1+EI.*qz2.*t4.*t34.*t67.*t68.*t81.*1.8e+1+EI.*qz3.*t4.*t34.*t67.*t69.*t81.*8.1e+1)-EA.*t241.*t324,qy1.*(t242+t252+t254+(EA.*qy1.*t2.*t32.*t70.*t323)./2.0+EA.*qy2.*t2.*t32.*t64.*t66.*t323+EA.*qy3.*t2.*t32.*t64.*t67.*t323.*(3.0./2.0)-EA.*qy1.*t3.*t23.*t34.*t64.*t65.*t241-EA.*qy2.*t3.*t23.*t34.*t65.*t66.*t241-EA.*qy2.*t3.*t23.*t34.*t64.*t68.*t241.*2.0-EA.*qy3.*t3.*t23.*t34.*t65.*t67.*t241.*(3.0./2.0)-EA.*qy3.*t3.*t23.*t34.*t64.*t69.*t241.*(9.0./2.0))+qy2.*(t244+t251+t260+EA.*qy2.*t2.*t32.*t72.*t323.*2.0+EA.*qy1.*t2.*t32.*t64.*t66.*t323+EA.*qy3.*t2.*t32.*t66.*t67.*t323.*3.0-EA.*qy1.*t3.*t23.*t34.*t65.*t66.*t241-EA.*qy1.*t3.*t23.*t34.*t64.*t68.*t241.*2.0-EA.*qy2.*t3.*t23.*t34.*t66.*t68.*t241.*8.0-EA.*qy3.*t3.*t23.*t34.*t66.*t69.*t241.*9.0-EA.*qy3.*t3.*t23.*t34.*t67.*t68.*t241.*6.0)+qy3.*(t246+t253+t259+EA.*qy3.*t2.*t32.*t73.*t323.*(9.0./2.0)+EA.*qy1.*t2.*t32.*t64.*t67.*t323.*(3.0./2.0)+EA.*qy2.*t2.*t32.*t66.*t67.*t323.*3.0-EA.*qy1.*t3.*t23.*t34.*t65.*t67.*t241.*(3.0./2.0)-EA.*qy1.*t3.*t23.*t34.*t64.*t69.*t241.*(9.0./2.0)-EA.*qy2.*t3.*t23.*t34.*t66.*t69.*t241.*9.0-EA.*qy2.*t3.*t23.*t34.*t67.*t68.*t241.*6.0-EA.*qy3.*t3.*t23.*t34.*t67.*t69.*t241.*2.7e+1)+qz1.*(t243+t256+t258+(EA.*qz1.*t2.*t32.*t70.*t323)./2.0+EA.*qz2.*t2.*t32.*t64.*t66.*t323+EA.*qz3.*t2.*t32.*t64.*t67.*t323.*(3.0./2.0)-EA.*qz1.*t3.*t23.*t34.*t64.*t65.*t241-EA.*qz2.*t3.*t23.*t34.*t65.*t66.*t241-EA.*qz2.*t3.*t23.*t34.*t64.*t68.*t241.*2.0-EA.*qz3.*t3.*t23.*t34.*t65.*t67.*t241.*(3.0./2.0)-EA.*qz3.*t3.*t23.*t34.*t64.*t69.*t241.*(9.0./2.0))+qz2.*(t245+t255+t262+EA.*qz2.*t2.*t32.*t72.*t323.*2.0+EA.*qz1.*t2.*t32.*t64.*t66.*t323+EA.*qz3.*t2.*t32.*t66.*t67.*t323.*3.0-EA.*qz1.*t3.*t23.*t34.*t65.*t66.*t241-EA.*qz1.*t3.*t23.*t34.*t64.*t68.*t241.*2.0-EA.*qz2.*t3.*t23.*t34.*t66.*t68.*t241.*8.0-EA.*qz3.*t3.*t23.*t34.*t66.*t69.*t241.*9.0-EA.*qz3.*t3.*t23.*t34.*t67.*t68.*t241.*6.0)+qz3.*(t247+t257+t261+EA.*qz3.*t2.*t32.*t73.*t323.*(9.0./2.0)+EA.*qz1.*t2.*t32.*t64.*t67.*t323.*(3.0./2.0)+EA.*qz2.*t2.*t32.*t66.*t67.*t323.*3.0-EA.*qz1.*t3.*t23.*t34.*t65.*t67.*t241.*(3.0./2.0)-EA.*qz1.*t3.*t23.*t34.*t64.*t69.*t241.*(9.0./2.0)-EA.*qz2.*t3.*t23.*t34.*t66.*t69.*t241.*9.0-EA.*qz2.*t3.*t23.*t34.*t67.*t68.*t241.*6.0-EA.*qz3.*t3.*t23.*t34.*t67.*t69.*t241.*2.7e+1)+qy1.*(t82+t89+t91+EI.*qy1.*t5.*t23.*t36.*t64.*t65+EI.*qy2.*t5.*t23.*t36.*t65.*t66.*4.0+EI.*qy2.*t5.*t23.*t36.*t64.*t68.*2.0+EI.*qy3.*t5.*t23.*t36.*t65.*t67.*(2.7e+1./2.0)+EI.*qy3.*t5.*t23.*t36.*t64.*t69.*(9.0./2.0))+qy2.*(t84+t88+t97+EI.*qy1.*t5.*t23.*t36.*t65.*t66.*4.0+EI.*qy1.*t5.*t23.*t36.*t64.*t68.*2.0+EI.*qy2.*t5.*t23.*t36.*t66.*t68.*3.2e+1+EI.*qy3.*t5.*t23.*t36.*t66.*t69.*3.6e+1+EI.*qy3.*t5.*t23.*t36.*t67.*t68.*5.4e+1)+qy3.*(t86+t90+t96+EI.*qy1.*t5.*t23.*t36.*t65.*t67.*(2.7e+1./2.0)+EI.*qy1.*t5.*t23.*t36.*t64.*t69.*(9.0./2.0)+EI.*qy2.*t5.*t23.*t36.*t66.*t69.*3.6e+1+EI.*qy2.*t5.*t23.*t36.*t67.*t68.*5.4e+1+EI.*qy3.*t5.*t23.*t36.*t67.*t69.*2.43e+2)+qz1.*(t83+t93+t95+EI.*qz1.*t5.*t23.*t36.*t64.*t65+EI.*qz2.*t5.*t23.*t36.*t65.*t66.*4.0+EI.*qz2.*t5.*t23.*t36.*t64.*t68.*2.0+EI.*qz3.*t5.*t23.*t36.*t65.*t67.*(2.7e+1./2.0)+EI.*qz3.*t5.*t23.*t36.*t64.*t69.*(9.0./2.0))+qz2.*(t85+t92+t99+EI.*qz1.*t5.*t23.*t36.*t65.*t66.*4.0+EI.*qz1.*t5.*t23.*t36.*t64.*t68.*2.0+EI.*qz2.*t5.*t23.*t36.*t66.*t68.*3.2e+1+EI.*qz3.*t5.*t23.*t36.*t66.*t69.*3.6e+1+EI.*qz3.*t5.*t23.*t36.*t67.*t68.*5.4e+1)+qz3.*(t87+t94+t98+EI.*qz1.*t5.*t23.*t36.*t65.*t67.*(2.7e+1./2.0)+EI.*qz1.*t5.*t23.*t36.*t64.*t69.*(9.0./2.0)+EI.*qz2.*t5.*t23.*t36.*t66.*t69.*3.6e+1+EI.*qz2.*t5.*t23.*t36.*t67.*t68.*5.4e+1+EI.*qz3.*t5.*t23.*t36.*t67.*t69.*2.43e+2)+EA.*t241.*t323];