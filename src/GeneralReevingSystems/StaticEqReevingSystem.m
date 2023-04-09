function Eq = StaticEqReevingSystem(qlam,reevSys)

t = 0.0;

n = reevSys.n;

%Extract independent coordinates and velocities
q = qlam(1:n,1);
lam = qlam(n+1:end,1);

v = zeros(n,1);

p = CalculateAllCoordinates(t,q,reevSys);
Cp = JacobianConstraints(p,reevSys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATION OF DEPENDENT VELOCITIES IN v %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = matrixB(q,reevSys);
D = matrixD(t,reevSys);

if isempty(Cp)
    Cq = [];
else
    Cq = Cp*B;
end

dp = B*v + D;

%%%%%%%%%%%%%%%%
%%% DYNAMICS %%%
%%%%%%%%%%%%%%%%

% Calculation of equations of motion

[M_p,Q_p] = MassMatrixAppliedForcesReevingSystem(p,dp,reevSys);

Eq = zeros(length(qlam),1);

Qap = B'*Q_p;

if isempty(Cp)
    Eq(1:n,1) = Qap;
else
    Eq(1:n,1) = Cq'*lam-Qap;
    Eq(n+1:end,1) = Constraints(p,reevSys);
end


