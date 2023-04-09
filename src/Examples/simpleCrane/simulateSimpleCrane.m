restoredefaultpath;
addpath("../../GeneralReevingSystems");
addpath("../../WireRopeForces");
addpath("../../mexInterface");
addpath("RB_kinematics");

% Getting system parameters %
reevSys = CraneParams;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinates: [r2x r2y r2z tet2_0 tet2_1 tet2_2 tet2_3 sa2 qbx1 qbx2 qbx3 qby1 qbz1] %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference value of coordinates (undeformed configuration) %
qref = [reevSys.d+reevSys.R 0.0 -(reevSys.l0 + reevSys.b) 1.0 0.0 0.0 0.0 reevSys.d 0.0 0.0 0.0 0.0 0.0]';
n = reevSys.n;

% Estimation of static displacement %
delta = reevSys.m2*9.81/(reevSys.EA/(reevSys.l0+reevSys.d));  % Elongation of both ropes
qref(3) = qref(3) - delta;  

% Initial guess of Lagrange multipliers associated with axial load
% constraint and Euler parameters %
lam = [0 0]';

% Calculation of static equilibrium position %
if strcmp(reevSys.NLeq_method,'NewtonRaphson')
    q0lam = NewtonRaphson(@StaticEqReevingSystem,[qref' lam']',reevSys);
elseif strcmp(reevSys.NLeq_method,'fsolve')
    q0lam = fsolve(@(qlam) StaticEqReevingSystem(qlam,reevSys),[qref' lam']',optimset('TolFun',1e-3));
end

% Generalized coordinates in static equilibrium %
qst = q0lam(1:n,1);
% Initial value of independent coordinates %
qind0 = q0lam(reevSys.Ind);
% Initial guess of dependent coordinates %
reevSys.qdep0 = q0lam(reevSys.Dep);

% Calculation of all coordinates in static equilibrium %
pst = CalculateAllCoordinates(0.0,qst,reevSys);

% Calculation of axial force field along rope %
sa = 0:0.01:qst(8);
for i = 1:length(sa)
    Pa(i) = AxialLoad(pst(7+1:7+17),sa(i),reevSys.EA_wr(1),reevSys);
end
sb = qst(8):0.01:(reevSys.d + reevSys.l0);
for i = 1:length(sb)
    Pb(i) = AxialLoad(pst(7+17+1:7+17+17),sb(i),reevSys.EA_wr(2),reevSys);
end

% Plots axial force field %
figure(10);
plot(sa,0.001*Pa,'b',sb,0.001*Pb,'g');
title('Axial load on wire ropes. Initial static position')
xlabel('Arc-length along rope (m)');
ylabel('Axial load (KN)');
legend('Horizontal element a','Vertical element b');

%%%%%%%% 
% Simulation %
%%%%%%%%

% Initial value of independent velocities %
vind0 = zeros(size(reevSys.Ind));

if strcmp(reevSys.ODE_method,'ode15s')
    tspan = 0:0.001:5.0;
    options = odeset('MaxStep',0.001);
    [t y] = ode15s(@(t,y) EqMotReevingSystem(t,y,reevSys),tspan, [qind0' vind0']',options);
elseif strcmp(reevSys.ODE_method,'ode45')
     tspan = 0:0.001:5.0;
    options = odeset('MaxStep',0.001);
    [t y] = ode45(@(t,y) EqMotReevingSystem(t,y,reevSys),tspan, [qind0' vind0']',options);
elseif strcmp(reevSys.ODE_method,'RungeKutta4')
    [t y] = RungeKutta4(@EqMotReevingSystem,[0 5.0], [qind0' vind0']', 0.001, reevSys);
end

%%%%%%%%%%%%%%%%%
%%% POSTPROCESSING %%%
%%%%%%%%%%%%%%%%%

% Calculation of all coordinates all instants %
qdep = reevSys.qdep0;

p = zeros(reevSys.np,length(t));
% Axial load at end '1' of element 'a' %
P1 = zeros(length(t),1);
% Axial load at end '2' of element 'b' %
P2 = zeros(length(t),1);

for i = 1:length(t)
    [p(:,i), Cp, qdep] = Calculate_p_from_qind(t(i),y(i,1:length(reevSys.Ind))',qdep,reevSys);
    P1(i) = AxialLoad(p(7+1:7+17,i),p(7+16,i),reevSys.EA_wr(1),reevSys); 
    P2(i) = AxialLoad(p(7+17+1:7+17+17,i),p(7+17+17,i),reevSys.EA_wr(2),reevSys); 
end

% Plot time-history of axial loads %
figure(20);
plot(t,0.001*P1,'b',t,0.001*P2,'g')
title('Axial load on wire ropes')
xlabel('Time (s)');
ylabel('Axial load (KN)');
legend('End 1 of element a','End 2 of element b');

% Animation of crane motion %
% Number of frames of the animation %
N = 100;

% Resamples the data %
t2 = t(1):(t(end)-t(1))/(N-1):t(end);
p2 = zeros(reevSys.np,N);
for i = 1:reevSys.np
    p2(i,:) = interp1(t,p(i,:),t2);
end

% Animates the motion of the system %
for i = 1:length(t2)
    Fot(i) = AnimateCrane(t2(i),p2(:,i),reevSys);
end

rmpath("../../GeneralReevingSystems");
rmpath("../../WireRopeForces");
rmpath("../../mexInterface");
rmpath("RB_kinematics");