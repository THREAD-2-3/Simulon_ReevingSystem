function P = getAxialLoads(p,reevSys)

P = zeros(reevSys.nwr,1);

ncwr = reevSys.ncwr;

for i = 1:reevSys.nwr
    qwr = p(reevSys.nrb+ncwr*(i-1)+1:reevSys.nrb+ncwr*i,1);
    P(i,1) = AxialLoad(qwr,qwr(ncwr),reevSys.EA_wr(i),reevSys); 
end