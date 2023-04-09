function C = ConstraintsOfqdep(qdep,reevSys)

t = reevSys.t;
qind = reevSys.qind;

q = zeros(reevSys.n,1);

q(reevSys.Ind,1) = qind;
q(reevSys.Dep,1) = qdep;

p = CalculateAllCoordinates(t,q,reevSys);
C = Constraints(p,reevSys);