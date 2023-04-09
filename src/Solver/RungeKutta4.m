function [t y] = RungeKutta4(fun, tspan, y0, At, reevSys)

t0 = tspan(1);
tfin = tspan(2);

t(1) = t0;
y(:,1) = y0;

i = 2;
tt = t0;
yy = y0;

while tt<=tfin
    [k1,P] = feval(fun,tt,yy, reevSys);
    k2 = feval(fun,tt+At/2,yy+(At/2)*k1, reevSys);
    k3 = feval(fun,tt+At/2,yy+(At/2)*k2, reevSys);
    k4 = feval(fun,tt+At,yy+At*k3, reevSys);
    
    yy = yy + (At/6)*(k1+2*k2+2*k3+k4);
    tt = tt+At;
    
    t(i) = tt;
    y(:,i) = yy;
    i = i+1;
  
    if(min(P)<0)
        [Pneg,j] = min(P);
        ropeName = {'la','lb','lc','ld','ra','rb','rc','rd'};
        disp('Negative tension found in rope:');
        disp(ropeName{j})
        disp('Simulation stopped');
        break
    end
end
    
    
    