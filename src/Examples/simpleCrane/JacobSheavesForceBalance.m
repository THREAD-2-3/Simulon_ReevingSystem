function Cp = JacobSheavesForceBalance(p,reevSys)

ncwr = reevSys.ncwr;

qa = p(7+1:7+ncwr,1);
qb = p(7+ncwr+1:7+2*ncwr,1);

Pqa = JacobAxialLoad(qa,qa(ncwr),reevSys.EA_wr(1),reevSys);
Pqb = JacobAxialLoad(qb,qb(ncwr-1),reevSys.EA_wr(2),reevSys);

Cp = [zeros(1,7) Pqa -Pqb];