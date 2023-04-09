function dv = VelocidadVector_3D(q,dq,v)

w = Omega(q,dq);
A = RotMat(q);

dv = cross(w,A*v);
