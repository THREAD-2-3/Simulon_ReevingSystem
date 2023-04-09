function L = getLengths(p,reevSys)

ncwr = reevSys.ncwr; % Number of coordinates per element

L = zeros(reevSys.nwr,1);

for i = 1:reevSys.nwr
    qwr = p(reevSys.nrb+ncwr*(i-1)+1:reevSys.nrb+ncwr*i,1);
    L(i,1) = sqrt((qwr(1)-qwr(4))^2+(qwr(2)-qwr(5))^2+(qwr(3)-qwr(6))^2);
end