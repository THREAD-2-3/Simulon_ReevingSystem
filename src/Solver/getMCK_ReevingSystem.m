function [M,C,K] = getMCK_ReevingSystem(t,q,v,reevSys)

n = reevSys.n;
K = zeros(n,n);
C = zeros(n,n);

p = CalculateAllCoordinates(t,q,reevSys);
B = matrixB(q,reevSys);
D = matrixD(t,reevSys);
dp = B*v + D;

M_p = massMatrixReevingSystem(p,reevSys);
M = B'*M_p*B;

incr = 1e-3;

for i = 1:n
    q(i) = q(i)+incr; p = CalculateAllCoordinates(t,q,reevSys); B = matrixB(q,reevSys);
    Qp_mas = appliedForcesReevingSystem(t,p,dp,reevSys);
    Q_mas = B'*Qp_mas;

    q(i) = q(i)-2*incr; p = CalculateAllCoordinates(t,q,reevSys); B = matrixB(q,reevSys);
    Qp_menos = appliedForcesReevingSystem(t,p,dp,reevSys);
    Q_menos = B'*Qp_menos;

    q(i) = q(i)+incr;
    K(:,i) = (1/(2*incr))*(Q_mas-Q_menos);
end

K = -K;

B = matrixB(q,reevSys);

for i = 1:n
    v(i) = v(i)+incr; dp = B*v + D;
    Qp_mas = appliedForcesReevingSystem(t,p,dp,reevSys);
    Q_mas = B'*Qp_mas;

    v(i) = v(i)-2*incr; dp = B*v + D;
    Qp_menos = appliedForcesReevingSystem(t,p,dp,reevSys);
    Q_menos = B'*Qp_menos;

    v(i) = v(i)+incr;
    C(:,i) = (1/(2*incr))*(Q_mas-Q_menos);
end

C = -C;