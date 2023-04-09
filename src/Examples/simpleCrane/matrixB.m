function B = matrixB(q,reevSys)

q2 = q(1:7,1);
A2 = rotMat(q2);
G2 = GLoc(q2);

n = reevSys.n;
nrb = reevSys.nrb;
ncwr = reevSys.ncwr;

B2 = zeros(nrb,n);
Ba = zeros(ncwr,n);
Bb = zeros(ncwr,n);

B2(1:7,1:7) = eye(7);

Ba(17,8) = 1.0;  % sa2

u_b2 = [0 0 reevSys.b]';
u_b2_sk = skew(u_b2);

Bb(4:6,1:3) = eye(3);
Bb(4:6,4:7) = -A2*u_b2_sk*G2;

Bb(7,9) = 1.0;   % qbx1
Bb(8,10) = 1.0;  % qbx2
Bb(9,11) = 1.0;  % qbx3
Bb(10,12) = 1.0;   % qby1
Bb(13,13) = 1.0;   % qbz1

Bb(16,8) = 1.0; % sb1

B = [B2' Ba' Bb']';