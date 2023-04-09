function Cq_RB = RigidBodyJacobianConstraints(p,reevSys)

% Euler parameters as unit quaternions

Cq_RB = zeros(1,length(p));

q2 = p(1:7,1);

Cq_RB(1,4) = 2*q2(4);
Cq_RB(1,5) = 2*q2(5);
Cq_RB(1,6) = 2*q2(6);
Cq_RB(1,7) = 2*q2(7);