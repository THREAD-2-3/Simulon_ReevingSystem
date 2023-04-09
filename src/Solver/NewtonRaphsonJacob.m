function [x, convergencia] = NewtonRaphsonJacob(fun,jac,x0,reevSys)

%PARAMETROS DEL NEWTONRAPHSON

% Newton-Raphson parameters

AbsErr=10^(-3);
MaxIter=100;

incr=10^(-5);

Niter=0;
ok=0;

while ((ok==0)&&(Niter<MaxIter))
    Niter=Niter+1;
    
    CC = feval(fun,x0,reevSys);
    Cq = feval(jac,x0,reevSys);
    
    % Evaluates convergency
    
    if (max(abs(CC))<AbsErr)
        ok=1;
    end
    
    % Newton-Raphson iteration step
    
    x=x0-Cq\CC;
    x0=x;
end
   
if Niter == MaxIter
    disp("Maximun number of iterations exceeded in Newton Raphson");
end



