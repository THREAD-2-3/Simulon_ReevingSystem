function [yp,P] = EqMotReevingSystem(t,y,reevSys)

t

nind = length(reevSys.Ind);
n = reevSys.n;
m = reevSys.m;

%Extract independent coordinates and velocities

qind = y(1:nind,1);
vind = y(nind+1:2*nind,1);

persistent qdep vdep t_before a_before  %Keeps the intial estimation of qdep

if(reevSys.m>0)
    if t == 0.0
        qdep = reevSys.qdep0;
        vdep = zeros(size(qdep));
        t_before = 0.0;
        a_before = zeros(n,1);
    end
    
    % Update initial guess
    
    qdep = qdep + vdep*(t-t_before);
    t_before = t;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CALCULATION OF DEPENDENT COORDINATES IN q %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    reevSys.a_before = a_before;
    [p, Cp, qdep] = Calculate_p_from_qind(t,qind,qdep,reevSys);
    
    q(reevSys.Ind,1) = qind;
    q(reevSys.Dep,1) = qdep;
    
else
    q(reevSys.Ind,1) = qind;
    p = CalculateAllCoordinates(t,q,reevSys);
    Cp = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATION OF DEPENDENT VELOCITIES IN v %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = matrixB(q,reevSys);
D = matrixD(t,reevSys);

% Jacobian of constraint with respect to q
if isempty(Cp)
    Cq = [];
    v = vind;
else
    Cq = Cp*B;
    Cqind = Cq(:,reevSys.Ind);
    Cqdep = Cq(:,reevSys.Dep);
    
    vdep = -Cqdep\(Cqind*vind);
    
    v(reevSys.Ind,1) = vind;
    v(reevSys.Dep,1) = vdep;
    
    Cqind = Cq(:,reevSys.Ind);
    Cqdep = Cq(:,reevSys.Dep);
    
    vdep = -Cqdep\(Cqind*vind);
    
    v(reevSys.Ind,1) = vind;
    v(reevSys.Dep,1) = vdep;
end

dp = B*v + D;

dB = dtmatrixB(q,v,reevSys);
dD = dtmatrixD(t,reevSys);

dCp = dtJacobianConstraints(p,dp,reevSys);

if isempty(dCp)
    dCq = [];
else
    dCq = dCp*B + Cp*dB;
end

%%%%%%%%%%%%%%%%
%%% DYNAMICS %%%
%%%%%%%%%%%%%%%%

% Calculation of equations of motion
[M_p,Q_p] = MassMatrixAppliedForcesReevingSystem(p,dp,reevSys);
% Adds motor torque (PI-control)
[Q_p, reevSys] = getMotorTorque(t,qind,vind,Q_p,reevSys);
% Adds sheave unbalance
[Q_p, reevSys] = getSheavesUnbalance(t,qind,vind,Q_p,reevSys);

M = B'*M_p*B;

if isempty(Cq)
    AA = zeros(n,n);
    AA(1:n,1:n) = M;
    
    bb = zeros(n,1);
    bb(1:n,1) = B'*(Q_p-M_p*(dB*v+dD));
else
    AA = zeros(n+m,n+m);
    AA(1:n,1:n) = M;
    AA(n+1:n+m,1:n) = Cq;
    AA(1:n,n+1:n+m) = Cq';
    
    bb = zeros(n+m,1);
    bb(1:n,1) = B'*(Q_p-M_p*(dB*v+dD));
    bb(n+1:n+m,1) = -dCq*v;
end
    
xx = AA\bb;
a = xx(1:n,1);

yp(1:nind,1) = v(reevSys.Ind);
yp(nind+1:2*nind,1) = a(reevSys.Ind);

a_before = a;

P = getAxialLoads(p,reevSys);

