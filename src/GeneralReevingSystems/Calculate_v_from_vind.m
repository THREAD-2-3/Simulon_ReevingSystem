function v = Calculate_v_from_vind(t,vind,q,Cp,reevSys)

if isempty(Cp)
    v = vind;
else
    B = matrixB(q,reevSys);
    D = matrixD(t,reevSys);
    
    % Jacobian of constraint with respect to q
    Cq = Cp*B;
    
    Cqind = Cq(:,reevSys.Ind);
    Cqdep = Cq(:,reevSys.Dep);
    
    vdep = -Cqdep\(Cqind*vind + Cp*D);
    
    v(reevSys.Ind,1) = vind;
    v(reevSys.Dep,1) = vdep;
end