function dr = VelocidadPunto_3D(r,q,dq)

w = Omega(q,dq);
A = RotMat(q);

dr = dq(1:3)+cross(w,A*r);
