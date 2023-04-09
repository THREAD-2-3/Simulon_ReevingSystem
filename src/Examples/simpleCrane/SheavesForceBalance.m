function C = SheavesForceBalance(p,reevSys)

ncwr = reevSys.ncwr;

qa = p(7+1:7+ncwr,1);
qb = p(7+ncwr+1:7+2*ncwr,1);

Pa = AxialLoad(qa,qa(ncwr),reevSys.EA_wr(1),reevSys);
Pb = AxialLoad(qb,qb(ncwr-1),reevSys.EA_wr(2),reevSys);

C = Pa-Pb;