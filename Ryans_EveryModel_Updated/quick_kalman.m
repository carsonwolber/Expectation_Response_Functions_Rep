%QUICK_KALMAN - Original representation of the Barsky-Sims (2012)
%information structure
function [out,GG0,GG1,GG2,F,sqQ,H,sqR,GG3] = quick_kalman(param)

rhoga = param.rhoga;
siges = param.siges;
sigga = param.sigga;
sigea = param.sigea;


F = [0,1;0 rhoga]; %state evol
H = [1 0; 0 1]';  %obs coeffs

Q = [sigea^2, 0; 0, sigga^2];      %shocks to growth
R = [0 0;0,siges^2]; %measurement errors

P = .0001*eye(2);
crit = 1;
while crit > 1e-15
    P_new = F*(P-P*H*((H'*P*H+R)\(H'*P)))*F' + Q;
    crit = max(max(abs(P-P_new)));
    P = P_new;
end

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
