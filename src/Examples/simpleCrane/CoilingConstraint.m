function [sa1,dsa1,ddsa1] = CoilingConstraint(t,reevSys)

if t < reevSys.t0
    sa1 = 0.5*reevSys.V*(t-(reevSys.t0/pi)*sin(pi*t/reevSys.t0));
    dsa1 = 0.5*reevSys.V*(1-cos(pi*t/reevSys.t0));
    ddsa1 = 0.5*reevSys.V*(pi/reevSys.t0)*(sin(pi*t/reevSys.t0));
else
    sa1 = 0.5*reevSys.V*reevSys.t0 + reevSys.V*(t-reevSys.t0);
    dsa1 = reevSys.V;
    ddsa1 = 0.0;
end
