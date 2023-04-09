function [M,Q] = MassMatrixAppliedForcesReevingSystem(p,dp,reevSys)

[M_RB, Q_RB] = MassMatrixAppliedForcesRigidBodies(p,dp,reevSys);
[M_WR, Q_WR] = MassMatrixAppliedForcesWireropeElements(p,dp,reevSys);

n = reevSys.nrb + reevSys.ncwr*reevSys.nwr + reevSys.nccr*reevSys.ncr;

M = zeros(n,n);
Q = zeros(n,1);

M(1:reevSys.nrb,1:reevSys.nrb) = M_RB;
M(reevSys.nrb+1:n,reevSys.nrb+1:n) = M_WR;

Q(1:reevSys.nrb,1) = Q_RB;
Q(reevSys.nrb+1:n,1) = Q_WR;

Qsusp = suspForces(p,dp,reevSys);

Q = Q + Qsusp;




