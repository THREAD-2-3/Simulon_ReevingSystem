function a = Calculate_a_from_aind(t,aind,v,p,Cp,reevSys)

if isempty(Cp)
    a = aind;
else
    q = p(reevSys.q_in_p);
    
    B = matrixB(q,reevSys);
    D = matrixD(t,reevSys);
    dB = dtmatrixB(q,v,reevSys);
    dD = dtmatrixD(t,reevSys);
    
    % Jacobian of constraint with respect to q
    Cq = Cp*B;
    
    % Time-derivative of Jacobian of constraint with respect to q
    dp = B*v+D;
    dCp = dtJacobianConstraints(p,dp,reevSys);
    dCq = dCp*B + Cp*dB;
    
    Cqind = Cq(:,reevSys.Ind);
    Cqdep = Cq(:,reevSys.Dep);
    
    adep = -Cqdep\(Cqind*aind + dCq*v + Cp*dD);
    
    a(reevSys.Ind,1) = aind;
    a(reevSys.Dep,1) = adep;   
end