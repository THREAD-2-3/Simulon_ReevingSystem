function D = matrixD(t,reevSys)

nrb = reevSys.nrb;
ncwr = reevSys.ncwr;

D2 = zeros(nrb,1);
Da = zeros(ncwr,1);
Db = zeros(ncwr,1);

% [sa1,dsa1,ddsa1] = CoilingConstraint(t,reevSys);
% [alpha,dalpha,ddalpha] = AzimutConstraint(t,reevSys);

dsa1 = interp1(reevSys.sa1_lookup(:,1),reevSys.sa1_lookup(:,3),t,'linear');
alpha = interp1(reevSys.alpha_lookup(:,1),reevSys.alpha_lookup(:,2),t,'linear');
dalpha = interp1(reevSys.alpha_lookup(:,1),reevSys.alpha_lookup(:,3),t,'linear');

d = reevSys.d;
R = reevSys.R;

Da(4:6,1) = [-d*sin(alpha)*dalpha d*cos(alpha)*dalpha 0]';

Da(ncwr-1,1) = dsa1;

Db(1:3,1) = [-(d+R)*sin(alpha)*dalpha (d+R)*cos(alpha)*dalpha 0]';

D = [D2' Da' Db']';

