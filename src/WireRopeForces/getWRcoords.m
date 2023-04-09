function qwr = getWRcoords(t,q,reevSys)

p = CalculateAllCoordinates(t,q,reevSys);
ncwr = reevSys.ncwr;
qwr = zeros(ncwr,reevSys.nwr);

for i = 1:reevSys.nwr
    qwr(:,i) = p(reevSys.nrb+ncwr*(i-1)+1:reevSys.nrb+ncwr*i,1);
end