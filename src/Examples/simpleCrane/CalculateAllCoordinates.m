function p = CalculateAllCoordinates(t,q,reevSys)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinates: [r2x r2y r2z tet2_0 tet2_1 tet2_2 tet2_3 sa2 qbx1 qbx2 qbx3] %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinates of rigid block

q2 = q(1:7,1);

% Coordinates of element 'a'

% [sa1,dsa1,ddsa1] = CoilingConstraint(t,reevSys);
% [alpha,dalpha,ddalpha] = AzimutConstraint(t,reevSys);

sa1 = interp1(reevSys.sa1_lookup(:,1),reevSys.sa1_lookup(:,2),t,'linear');
alpha = interp1(reevSys.alpha_lookup(:,1),reevSys.alpha_lookup(:,2),t,'linear');

ra1 = [0.0 0.0 reevSys.R]'; 
ra2 = [reevSys.d*cos(alpha) reevSys.d*sin(alpha) reevSys.R]'; 
qax = [0.0 0.0 0.0]';
qay = [0.0 0.0 0.0]';
qaz = [0.0 0.0 0.0]';

sa = [sa1 q(8)]';

qa = [ra1' ra2' qax' qay' qaz' sa']';

% Coordinates of element 'b'

A2 = rotMat(q2);
u_b2 = [0 0 reevSys.b]';

rb1 = [(reevSys.d+reevSys.R)*cos(alpha) (reevSys.d+reevSys.R)*sin(alpha) 0.0]'; 
rb2 = q2(1:3,1) + A2*u_b2;
qbx = [q(9) q(10) q(11)]'; 
qby = [q(12) 0.0 0.0]'; 
qbz = [q(13) 0.0 0.0]'; 

sb1 = q(8);
sb2 = reevSys.l0 + reevSys.d;

sb = [sb1 sb2]';

qb = [rb1' rb2' qbx' qby' qbz' sb']';

% Group all coordinates

p = [q2' qa' qb']';