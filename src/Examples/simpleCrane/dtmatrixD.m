function dD = dtmatrixD(t,reevSys)

nrb = reevSys.nrb;
ncwr = reevSys.ncwr;

dD2 = zeros(nrb,1);
dDa = zeros(ncwr,1);
dDb = zeros(ncwr,1);

% [sa1,dsa1,ddsa1] = CoilingConstraint(t,reevSys);
% [alpha,dalpha,ddalpha] = AzimutConstraint(t,reevSys);

ddsa1 = interp1(reevSys.sa1_lookup(:,1),reevSys.sa1_lookup(:,4),t,'linear');
alpha = interp1(reevSys.alpha_lookup(:,1),reevSys.alpha_lookup(:,2),t,'linear');
dalpha = interp1(reevSys.alpha_lookup(:,1),reevSys.alpha_lookup(:,3),t,'linear');
ddalpha = interp1(reevSys.alpha_lookup(:,1),reevSys.alpha_lookup(:,4),t,'linear');

d = reevSys.d;
R = reevSys.R;

dDa(4:6,1) = [-d*sin(alpha)*ddalpha-d*cos(alpha)*(dalpha^2) d*cos(alpha)*ddalpha-d*sin(alpha)*(dalpha^2) 0]';

dDa(ncwr-1,1) = ddsa1;

dDb(1:3,1) = [-(d+R)*sin(alpha)*ddalpha-(d+R)*cos(alpha)*(dalpha^2) (d+R)*cos(alpha)*ddalpha-(d+R)*sin(alpha)*(dalpha^2) 0]';

dD = [dD2' dDa' dDb']';






