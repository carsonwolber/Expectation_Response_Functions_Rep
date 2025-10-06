function [yerf,exitflag,GAM,ZZ,TQA,TQB]=gx_hx_erf(fy,fx,fyp,fxp,stake,yidx,xidx,xkidx,eqidx,nerf)

nx = length(xidx);
ny = length(yidx);

fyx = [fy,fx];
fyxp = [fyp,fxp];

A1 = fyxp(eqidx,xidx);
A2 = fyxp(eqidx,yidx);

B1 = -fyx(eqidx,xidx);
B2 = -fyx(eqidx,yidx);

Bx = -fyx(eqidx,xkidx);
Ax = -fyxp(eqidx,xkidx);

[S,T,Q,Z] = qz([A1,A2],[B1,B2]);
%stake = 1;
slt = (abs(diag(T))<stake*abs(diag(S)));  

nk=sum(slt);

%Reorder the system with stable eigs in upper-left
[S,T,Q,Z] = ordqz(S,T,Q,Z,slt);   


S11 = S(1:nx,1:nx);
S22 = S(nx+1:end,nx+1:end);
T22 = T(nx+1:end,nx+1:end);
Z11 = Z(1:nx,1:nx);
Z22 = Z(nx+1:end,nx+1:end);
Z12 = Z(1:nx,nx+1:end);
Z21 = Z(nx+1:end,1:nx);
Q2  = Q(nx+1:end,:);


%Catch cases with no/multiple solutions
exitflag = 1;
if nk>nx
    warning('The Equilibrium is Locally Indeterminate')
    exitflag=2;
elseif nk<nx
    warning('No Local Equilibrium Exists')
    exitflag = 0;
elseif rank(Z11)<nk
    warning('Invertibility condition violated')
    exitflag = 3;
end



yerf = zeros(ny,nerf);
T22_S22 = (T22\S22);
T22_S22_jj    = eye(size(T22));
T22_S22_jj_m1 = eye(size(T22));
T22_Q2_Bx = (T22\(Q2*Bx));
T22_Q2_Ax = (T22\(Q2*Ax));
ZZZ = (Z22-Z21*(Z11\Z12));
for jj = 0:nerf

    tmp = 0;
    if jj > 0
        %tmp = tmp - (Z22-Z21*(Z11\Z12))*(T22\S22)^(jj-1)*(T22\(Q2*Ax));  %the "ryan" term
        tmp = tmp - ZZZ*T22_S22_jj*T22_Q2_Ax;  %the "ryan" term
        T22_S22_jj = T22_S22_jj*T22_S22;
    end

    %tmp = -(Z22-Z21*(Z11\Z12))*(T22\S22)^(jj)*(T22\(Q2*Bx)); % the "kyle" term
    tmp = tmp-ZZZ*T22_S22_jj*T22_Q2_Bx ; % the "kyle" term

    yerf(:,jj+1) = real(tmp);
end

if nargout > 2
     GAM = T22\S22;
     TQB = (T22\(Q2*Bx));
     TQA = (T22\(Q2*Ax));
     ZZ  = (Z22-Z21*(Z11\Z12));
end