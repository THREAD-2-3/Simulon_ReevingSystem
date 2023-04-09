function dB = dtmatrixB(q,v,reevSys)

n = reevSys.n;
nrb = reevSys.nrb;
ncwr = reevSys.ncwr;

dB2 = zeros(nrb,n);
dBa = zeros(ncwr,n);
dBb = zeros(ncwr,n);

q2 = q(1:7,1);
v2 = v(1:7,1);

A2 = rotMat(q2);
G2 = GLoc(q2);
wL2 = G2*v2(4:7,1);
wL2_sk = skew(wL2);

u_b2 = [0 0 reevSys.b]';
u_b2_sk = skew(u_b2);

dBb(4:6,4:7) = -A2*wL2_sk*u_b2_sk*G2;

dB = [dB2' dBa' dBb']';






