function reevSys = CraneParams

% Wire rope elastic and inertia properties

r = 5e-3;               % Radius of wire rope cross-section
E = 2.1e11;             % Young modulus of wire steel 
A = pi*r^2;             % Area of cross-section    
I = (1/4)*pi*r^4;       % Moment of inertia fo cross-section 
roS = 7900;             % Volumetric density of steel  

reevSys.EA = E*A;
reevSys.EI = 0.3*E*I;     % Reduced by a factor 0.3 becuase the wire-rope cross-section does not behave as a solid cross-section 
reevSys.ro = roS*A;       % Mass per unit length
nu = 0.001;                     % Damping   

for i = 1:2
    reevSys.EA_wr(i) = reevSys.EA;
    reevSys.EI_wr(i) = 0.0;
    reevSys.ro_wr(i) = reevSys.ro;
    reevSys.nu_wr(i) = nu;
    reevSys.elementType{i} = 'wireRope';
end

% Geometric parameters of reeving system

reevSys.R = 0.1;
reevSys.l0 = 8;
reevSys.b = 0.3;
reevSys.d = 3;

% Rigid block inertia properties

reevSys.m2 = 1000;
reevSys.I2x = 2*reevSys.m2*(reevSys.b)^2/5;         % Like the moment of inertia of a sphere with R = b
reevSys.I2y = 2*reevSys.m2*(reevSys.b)^2/5;         % Like the moment of inertia of a sphere with R = b
reevSys.I2z = 2*reevSys.m2*(reevSys.b)^2/5;         % Like the moment of inertia of a sphere with R = b

% Wire rope coiling constraint

reevSys.V = 1;   
reevSys.t0 = 0.2;

% Crane rotation constraint

reevSys.Om = 0.3;  % CRANE WITH NO ROTATION
reevSys.t1 = 0.2;

% Motion lookup tables

t = 0:0.01:10;
t = t';

for i = 1:length(t)
     [alpha(i,1),dalpha(i,1),ddalpha(i,1)] = AzimutConstraint(t(i),reevSys);
     [sa1(i,1),dsa1(i,1),ddsa1(i,1)] = CoilingConstraint(t(i),reevSys);
end

% figure(445)
% plot(t,dalpha,'LineWidth',2);
% xlabel('Time (s)','FontSize', 14)
% ylabel('Velocity (m/s)','FontSize', 14)
% title('Payload target velocity','FontSize', 14)
% xlim([0 0.5])
% ylim([0 0.33])

reevSys.sa1_lookup = [t sa1 dsa1 ddsa1];
reevSys.alpha_lookup = [t alpha dalpha ddalpha];

% Position of independent and dependent coordinates in vector of coordinates q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinates: [r2x r2y r2z tet2_0 tet2_1 tet2_2 tet2_3 sa2 qbx1 qbx2 qbx3 qby1 qbz1] %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reevSys.n = 13;
reevSys.m = 2;
reevSys.g = reevSys.n-reevSys.m;   % Number of degrees of freedom

reevSys.Ind = [1:6 9:13]';       % Index of independent coordinates in q
reevSys.Dep = [7 8]';              % Index of dependent coordinates in q 

% Position of integration coordinates q in vector of coordinates p 

reevSys.q_in_p = [1:7 7+17 7+17+6+1 7+17+6+2 7+17+6+3 7+17+6+4 7+17+6+7]';

% Number of rigid body coordinates, number of wire-rope elements and number of coordinates per element

reevSys.nrb = 7;
reevSys.nwr = 2;
reevSys.ncwr = 17;
reevSys.ncr = 0;
reevSys.nccr = 14;

reevSys.np = reevSys.nrb + reevSys.nwr*reevSys.ncwr + reevSys.ncr*reevSys.nccr;

% Position of element coordinates in qWR

reevSys.elemCoords{1} = reevSys.nrb + (1:17);
reevSys.elemCoords{2} = reevSys.nrb + (17+1:2*17);

% Position of Rigid Body coordinates and Wire Rope coordinates in p 

reevSys.RBcoord = [1:7]';
reevSys.WRcoord = [8:41]';

% Number of modal coordinates

reevSys.mnX(1) = 0; % Number of modal coordinates of element 1 in X direction...
reevSys.mnX(2) = 3;
reevSys.mnY(1) = 0; % Number of modal coordinates of element 1 in Y direction...
reevSys.mnY(2) = 1;
reevSys.mnZ(1) = 0; 
reevSys.mnZ(2) = 1;

% Estimation of value of dependent coordinate

tet2_3_0 = 0;
sa20 = reevSys.d;
reevSys.qdep0 = [tet2_3_0 sa20]';

% Method to solve nonlinear algebraic equations
reevSys.NLeq_method = 'fsolve';
%reevSys.NLeq_method = 'NewtonRaphson';
%reevSys.NLeq_method = 'LMFsolve';

% Method to solve ODE
reevSys.ODE_method = 'ode15s';
% reevSys.ODE_method = 'ode45';
% reevSys.ODE_method = 'genAlpha';
% reevSys.ODE_method = 'RungeKutta4';