function x = NewtonRaphson(fun,x0,reevSys)

% Newton-Raphson parameters

AbsErr = 10^(-3);
MaxIter = 200;

incr=10^(-5);

Niter=0;
ok=0;

while ((ok==0)&&(Niter<MaxIter))
    
    Niter=Niter+1;
   
% Numerical calculation of Jacobian matrix    
    
    CC = feval(fun,x0,reevSys);
    
    for i=1:length(x0)  
        
        x0(i)=x0(i)+incr;
        CCi=feval(fun,x0,reevSys);
        
        Cq(:,i)=(CCi-CC)/incr;
        x0(i)=x0(i)-incr;
    end    
    
% Evaluates convergency  

    if (max(abs(CC))<AbsErr)
        ok=1;
    end     
    
%     CC
%     pause
%     Cq 
%     size(Cq)
%     rank(Cq)
%     pause
    
% Newton-Raphson iteration step  
    x=x0-Cq\CC;
    x0=x;
end


if Niter == MaxIter
    disp("Maximun number of iterations exceeded in Newton Raphson");
end