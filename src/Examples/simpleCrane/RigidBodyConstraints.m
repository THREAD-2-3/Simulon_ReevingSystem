function C_RB = RigidBodyConstraints(p,reevSys)

q2 = p(1:7,1);

% Euler parameters as unit quaternions

C_RB = q2(4)^2+q2(5)^2+q2(6)^2+q2(7)^2-1;