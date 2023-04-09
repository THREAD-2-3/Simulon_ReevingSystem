function dCq_RB = RigidBodydtJacobianConstraints(p,dp,reevSys)

dCq_RB = zeros(1,length(p));

dq2 = dp(1:7,1);

% Euler parameters as unit quaternions

dCq_RB(1,4) = 2*dq2(4);
dCq_RB(1,5) = 2*dq2(5);
dCq_RB(1,6) = 2*dq2(6);
dCq_RB(1,7) = 2*dq2(7);