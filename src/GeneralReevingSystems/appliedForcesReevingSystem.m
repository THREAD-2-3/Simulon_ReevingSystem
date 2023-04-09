function Q_p = appliedForcesReevingSystem(t,p,dp,reevSys)

% Calculation of equations of motion
[M_p,Q_p] = MassMatrixAppliedForcesReevingSystem(p,dp,reevSys);

q = p(reevSys.q_in_p);
dq = dp(reevSys.q_in_p);
qind = q(reevSys.Ind);
dqind = dq(reevSys.Ind);

% Adds motor torque (PI-control)
[Q_p, reevSys] = getMotorTorque(t,qind,dqind,Q_p,reevSys);
% Adds sheave unbalance
[Q_p, reevSys] = getSheavesUnbalance(t,qind,dqind,Q_p,reevSys);
