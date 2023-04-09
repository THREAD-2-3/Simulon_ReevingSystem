function dCp = dtJacobSheavesForceBalance(p,v,reevSys)

ncwr = reevSys.ncwr;

qa = p(7+1:7+ncwr,1);
qb = p(7+ncwr+1:7+2*ncwr,1);

va = v(7+1:7+ncwr,1);
vb = v(7+ncwr+1:7+2*ncwr,1);

dPqa = DtJacobAxialLoad(qa,va,qa(ncwr),reevSys.EA_wr(1),reevSys);
dPqb = DtJacobAxialLoad(qb,vb,qb(ncwr-1),reevSys.EA_wr(1),reevSys);

dCp = [zeros(1,7) dPqa -dPqb];