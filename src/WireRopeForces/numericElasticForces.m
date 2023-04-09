function Qel = numericElasticForces(qe,dqe,EA,EI,nu,reevSys)

ncwr = reevSys.ncwr;

s1 = qe(ncwr-1); s2 = qe(ncwr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integral of sextic polynomials with Gauss quadrature %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% shi1 = -0.932469514203152; shi2 = -0.661209386466265; shi3 = -0.238619186083197;
% w1 = 0.171324492379170; w2 = 0.360761573048139; w3 = 0.467913934572691;
% shi4 = -shi3; shi5 = -shi2; shi6 = -shi1;
% w4 = w3; w5 = w2; w6 = w1;
% 
% s_1 = 0.5*(s2-s1)*shi1+0.5*(s2+s1); 
% s_2 = 0.5*(s2-s1)*shi2+0.5*(s2+s1); 
% s_3 = 0.5*(s2-s1)*shi3+0.5*(s2+s1); 
% s_4 = 0.5*(s2-s1)*shi4+0.5*(s2+s1); 
% s_5 = 0.5*(s2-s1)*shi5+0.5*(s2+s1); 
% s_6 = 0.5*(s2-s1)*shi6+0.5*(s2+s1); 
% 
% d_udef1 = jacobDiffDefEnergy(qe,s_1,EA,EI);
% d_udef2 = jacobDiffDefEnergy(qe,s_2,EA,EI);
% d_udef3 = jacobDiffDefEnergy(qe,s_3,EA,EI);
% d_udef4 = jacobDiffDefEnergy(qe,s_4,EA,EI);
% d_udef5 = jacobDiffDefEnergy(qe,s_5,EA,EI);
% d_udef6 = jacobDiffDefEnergy(qe,s_6,EA,EI);
% 
% dUdef = 0.5*(s2-s1)*(w1*d_udef1+w2*d_udef2+w3*d_udef3+w4*d_udef4+w5*d_udef5+w6*d_udef6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeated midpoint formula %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 4;
h = (s2-s1)/n;

d_udef = zeros(ncwr,n+1);

for i = 1:n+1
    s = s1 + (i-1)*h;
    if nu == 0
        d_udef(:,i) = jacobDiffDefEnergy_noDamp(qe,s,EA,EI);
    else
        d_udef(:,i) = jacobDiffDefEnergy_Damp(qe,dqe,s,EA,EI,nu);
    end
end

dUdef = (h/2)*(d_udef(:,1)+d_udef(:,n+1));

for i = 2:n
    dUdef = dUdef + h*d_udef(:,i);
end

%%%%%%%%%%%%%%%%%
% Leibniz terms %
%%%%%%%%%%%%%%%%%

% These terms have to be included no matter what quadrature you use above %

e1 = zeros(size(qe));
e2 = zeros(size(qe));
e1(end-1,1) = 1;
e2(end,1) = 1;

if nu == 0
    dUdef = dUdef + diffDefEnergy_noDamp(qe,s2,EA,EI)*e2 - diffDefEnergy_noDamp(qe,s1,EA,EI)*e1;
else
    dUdef = dUdef + diffDefEnergy_Damp(qe,dqe,s2,EA,EI,nu)*e2 - diffDefEnergy_Damp(qe,dqe,s1,EA,EI,nu)*e1;
end

Qel = -dUdef;



