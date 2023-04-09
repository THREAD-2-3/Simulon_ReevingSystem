function [alpha,dalpha,ddalpha] = AzimutConstraint(t,reevSys)

if t < reevSys.t1
    alpha = 0.5*reevSys.Om*(t-(reevSys.t1/pi)*sin(pi*t/reevSys.t1));
    dalpha = 0.5*reevSys.Om*(1-cos(pi*t/reevSys.t1));
    ddalpha = 0.5*reevSys.Om*(pi/reevSys.t1)*(sin(pi*t/reevSys.t1));
else
    alpha = 0.5*reevSys.Om*reevSys.t1 + reevSys.Om*(t-reevSys.t1);
    dalpha = reevSys.Om;
    ddalpha = 0.0;
end
