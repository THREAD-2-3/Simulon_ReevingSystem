function [p, Cp, qdep] = Calculate_p_from_qind(t,qind,qdep,reevSys)

if strcmp(reevSys.NLeq_method,'NewtonRaphson')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matlab 'fsolve' to calculate dependent coordines %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    reevSys.qind = qind;
    reevSys.t = t;
    
    if reevSys.m ~= 0
        qdep = NewtonRaphson(@ConstraintsOfqdep,qdep,reevSys);
    end

    q = zeros(reevSys.n,1);
    q(reevSys.Ind,1) = qind;
    q(reevSys.Dep,1) = qdep;
    p = CalculateAllCoordinates(t,q,reevSys);
    Cp = JacobianConstraints(p,reevSys);
    
end

if strcmp(reevSys.NLeq_method,'LevenbergMarquardtFletcher')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matlab 'fsolve' to calculate dependent coordines %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    reevSys.qind = qind;
    reevSys.t = t;
    
    if reevSys.m ~= 0
        f = @(qdep) ConstraintsOfqdep(qdep,reevSys); % function of dummy variable x
        qdep = LMFsolve(f,qdep);
    end
    
    q = zeros(reevSys.n,1);
    q(reevSys.Ind,1) = qind;
    q(reevSys.Dep,1) = qdep;
    
    p = CalculateAllCoordinates(t,q,reevSys);
    Cp = JacobianConstraints(p,reevSys);
    
end

if strcmp(reevSys.NLeq_method,'fsolve')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matlab 'fsolve' to calculate dependent coordines %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    reevSys.qind = qind;
    reevSys.t = t;
    
    if reevSys.m ~= 0
        f = @(qdep) ConstraintsOfqdep(qdep,reevSys); % function of dummy variable x
        qdep = fsolve(f,qdep,optimset('MaxFunEvals',100,'MaxIter',100,'TolFun',1e-3,'Display','off'));
    end
    
    q = zeros(reevSys.n,1);
    q(reevSys.Ind,1) = qind;
    q(reevSys.Dep,1) = qdep;
    
    p = CalculateAllCoordinates(t,q,reevSys);
    Cp = JacobianConstraints(p,reevSys);
end

