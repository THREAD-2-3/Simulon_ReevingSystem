function P = getAxialLoadsWithModes(p,reevSys)

nn = 3;

P = zeros(reevSys.nwr,nn);

ncwr = reevSys.ncwr;

for i = 1:reevSys.nwr
    qwr = p(reevSys.nrb+ncwr*(i-1)+1:reevSys.nrb+ncwr*i,1);
    s = linspace(qwr(ncwr-1,1),qwr(ncwr,1),nn);
    for j = 1:nn
        P(i,j) = AxialLoad(qwr,s(j),reevSys.EA_wr(i),reevSys);
    end
%     figure(i);
%     plot(s,P(i,:));
end