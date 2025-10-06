%QUICK_KALMAN_NOISE - Noise representation of the Barsky-Sims (2012)
%information structure

function [out,H,F,Q,GG0,GG1,GG2,sqQ,P] = quick_kalman_noise(param)

rhoga = param.rhoga;
siges = param.siges;
sigga = param.sigga;
sigea = param.sigea;

%From Kyle's notes
delt = ((1+ rhoga^2 + sigga^2/sigea^2) - sqrt((1+rhoga^2+sigga^2/sigea^2)^2-4*rhoga^2))/(2*rhoga);
bet =   ((1+rhoga^2+sigga^2/sigea^2*(sigea^2+siges^2)/siges^2) - sqrt((1+rhoga^2+sigga^2/sigea^2*(sigea^2+siges^2)/siges^2)^2-4*rhoga^2))/(2*rhoga);


phi1 = rhoga+delt;
phi2 = -rhoga*delt;

lam0 = -rhoga*(sigea^2/sigga^2);
lam1 = -(1+delt^2)/delt*lam0;
lam2 = lam0;

sigex = sqrt(delt/rhoga)*sigga^2/sigea;  %double check sigea power
sigev = sqrt(delt/bet)*siges;



F = [phi1 phi2 0 0 0; 1 0 0 0 0; 0 1 0 0 0; 0 0 0  delt -bet; 0 0 0 0 0]; %state evol
H = [lam0 lam1 lam2 0 0; 1 0 0 1 0]';  %obs coeffs

Q = [sigex^2 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 sigev^2 sigev^2; 0 0 0 sigev^2 sigev^2];      %shocks to stats
R = zeros(2); %measurement errors

P = .0001*eye(5);
crit = 1;
tt = 1;
while crit > 1e-16 && tt<10000;
    P_new = F*(P-P*H*((H'*P*H+R)\(H'*P)))*F' + Q;
    crit = max(max(abs(P_new-P)));
    P = P_new;
    tt = tt+1;
end

%rcond(P)

%Depdendence on xi(t-1)
GG0 = P*H*((H'*P*H+R)\(H'*F));

%Depedence on xi(t-1,t-1)
GG1 = F + P*H*((H'*P*H+R)\(-H'*F));

%Depepdence on v(t), state shocks
GG2 = P*H*((H'*P*H+R)\(H'))*sqrt(Q);

%Depepdence on w(t), obs shocks
GG3 = P*H/(H'*P*H+R)*sqrt(R);

sqQ = sqrt(Q);
sqR = sqrt(R);

out = [GG0(:);GG1(:);GG2(:);GG3(:);F(:);sqQ(:)];
