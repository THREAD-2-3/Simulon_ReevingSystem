function [M,Q] = MassMatrixAppliedForcesRigidBodies(q,v,reevSys)

M = eye(7);

M(1:3,1:3) = reevSys.m2*eye(3);

q2 = q(1:7,1);
v2 = v(1:7,1);
A2 = rotMat(q2);
G2 = GLoc(q2);
wL2 = G2*v2(4:7,1);
wL2_sk = skew(wL2); 

I = eye(3);
I(1,1) = reevSys.I2x;
I(2,2) = reevSys.I2y;
I(3,3) = reevSys.I2z;

M(4:7,4:7) = G2'*I*G2; 

Qgrav = [0 0 -reevSys.m2*9.81 0 0 0 0]';

Qv = zeros(7,1);
Qv(4:7,1) = -G2'*wL2_sk*I*wL2;

Q = Qv+Qgrav;