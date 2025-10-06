function [Vy,mean_horizon] = erf_filter(GAM,ZZ,TQA,TQB,lam)
% -------------------------------------------------------------------------
% Compute variance decomposition over frequencies lam = [lam1,lam2] when
% 
% ytil = -ZZ*(I-GAM*L)^-1*(TQB+TQA*L)*eps(t)
%
% -------------------------------------------------------------------------

% Options: default 6-32 periods
if nargin < 5
    lam = [2*pi/32,2*pi/6];
end 


% Compute filtered variances due to each shock
lam1  = lam(1);
lam2  = lam(2);
II    = eye(size(GAM));


fy = @(w) spec_dens(w,GAM,ZZ,TQA,TQB,II);
Vy = 2.*real(integral(fy,lam1,lam2,'arrayvalued',true))'    ;


function out = spec_dens(w,GAM,ZZ,TQA,TQB,II)

T1 = -ZZ*((II-GAM.*exp(-1i.*w))\(TQB + TQA.*exp(-1i.*w)));
out = 1/(2*pi)*(T1*T1');
