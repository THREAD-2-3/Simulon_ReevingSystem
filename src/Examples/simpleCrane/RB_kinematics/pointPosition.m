function r = PosicionPunto_3D(r1,q)

A = RotMat(q);

r = q(1:3)+A*r1;