function [t,q,dq,ddq,lam] = generalizedAlpha(tspan, q0, dq0, reevSys)

t0 = tspan(1);
tfin = tspan(end);
At = tspan(2)-tspan(1);
N = length(t0:At:tfin);

n = reevSys.n;
m = reevSys.m;

alf_m = reevSys.alf_m;
alf_f = reevSys.alf_f;
bet = reevSys.bet;
gam = reevSys.gam;
tol = reevSys.tol;

betp = (1-alf_m)/(At^2*bet*(1-alf_f));
gamp = gam/(At*bet);

t = zeros(1,N);
q = zeros(n,N);
dq = zeros(n,N);
ddq = zeros(n,N);
lam = zeros(m,N);

ddq0 = zeros(n,1);
a = zeros(n,1);

t(1) = t0;
q(:,1) = q0;
dq(:,1) = dq0;

M = reevSys.M_st;
C = reevSys.C_st;
K = reevSys.K_st;

St = zeros(n+m,n+m);
bb = zeros(n+m,1);

DL = eye(n+m,n+m);
DR = eye(n+m,n+m);

DL(1:n,1:n) = At^2*bet*DL(1:n,1:n);
DR(n+1:n+m,n+1:n+m) = (1/(At^2*bet))*DR(n+1:n+m,n+1:n+m);

St(1:n,1:n) = betp*M + gamp*C + K; 

for i = 2:N

    t0 = t0+At

    % SOLVE EOM %

    q1 = q0 + At*dq0 + At^2*(0.5-bet)*a;
    dq1 = dq0 + At*(1-gam)*a;
    lam1 = zeros(m,1);

    a = (1/(1-alf_m))*(alf_f*ddq0 - alf_m*a);

    q1 = q1 + At^2*bet*a;
    dq1 = dq1 + At*gam*a;
    ddq1 = zeros(n,1);

    % Keeps values in case jacobians have to be updated %
    q10 = q1;
    dq10 = dq1;
    ddq10 = zeros(n,1);
    lam10 = zeros(m,1);

    for j = 1:19

        p = CalculateAllCoordinates(t0,q1,reevSys);
        B = matrixB(q1,reevSys);
        D = matrixD(t0,reevSys);
        dp = B*dq1 + D;

        M_p = massMatrixReevingSystem(p,reevSys);
        M1 = B'*M_p*B;
        Q_p = appliedForcesReevingSystem(t0,p,dp,reevSys);
        Q1 = B'*Q_p;

        if reevSys.m > 0
            C1 = Constraints(p,reevSys);
            Cp = JacobianConstraints(p,reevSys);
            Cq1 = Cp*B;

            r_lam = C1;
            r_q = M1*ddq1+Cq1'*lam1-Q1;
        else
            r_lam = [];
            C1 = [];
            Cq1 = [];
             r_q = M1*ddq1-Q1;
        end

        r = sqrt(norm(r_q)^2+norm(r_lam)^2);
        if r<tol
            break
        end
 
        St(1:n,1:n) = betp*M + gamp*C + K;
        St(1:n,n+1:n+m) = Cq1';
        St(n+1:n+m,1:n) = Cq1;

        St_hat = DL*St*DR;
        r_hat = DL*[r_q' r_lam']';

% if j == 11
%     rank(St_hat)
%     St_hat
%     pause
% end

        x_hat = -St_hat\r_hat;

        Aq = x_hat(1:n,1);
        Alam = (1/(At^2*bet))*x_hat(n+1:n+m,1);

        q1 = q1 + Aq;
        dq1 = dq1 + gamp*Aq;
        ddq1 = ddq1 + betp*Aq;
        lam1 = lam1 + Alam;

        % Update jacobians %

        if ~mod(j,10)   % When j = 10 %
            q1 = q10;
            dq1 = dq10;
            ddq1 = ddq10;
            lam1 = lam10;

            [M,C,K] = getMCK_ReevingSystem(t0,q1,dq1,reevSys);
            disp('Re-calculation of Jacobians');
        end
    end

    if j == 19
        disp('Time step terminated without convergence');
    end

    % STORE VALUES AND UPDATE FOR NEXT STEP %

    t(i) = t0;
    q(:,i) = q1;
    dq(:,i) = dq1;
    ddq(:,i) = ddq1;
    lam(:,i) = lam1;

%     reevSys.fn = -lam1([3 5 7 9],1); % Update normal forces for the calculation of the tangential contact forces

    q0 = q1;
    dq0 = dq1;
    ddq0 = ddq1;

    a = a + ((1-alf_f)/(1-alf_m))*ddq1;
end