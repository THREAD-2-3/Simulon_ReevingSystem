function M_p = massMatrixReevingSystem(p,reevSys)

[M_p,Q_p] = MassMatrixAppliedForcesReevingSystem(p,zeros(size(p)),reevSys);

