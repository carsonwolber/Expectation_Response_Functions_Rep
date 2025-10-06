% MODEL_SS - Return the steady state of the model computed analytically
%
% usage:
% 
% [ss, parameters] =model_ss(param)


function [Yss,Xss,param,set] = jr_ss(param,set)


%Upack parameters object
param_unpack


%BEGIN_EXTRACT_HERE

rbar = 1/bet; %actual domestic interest rate
h = hbar;


k = h*((1/bet - 1 + delt0)/alph)^(1/(alph-1));
i = k*delt0;
gdp = k^alph*h^(1-alph);
d = dshr*gdp;
c = gdp - i - d*(1-1/rbar);
x = c;
psii = (1-alph)*k^alph*h^(1-alph)/(h^thet*(c*thet-(1-alph)*k^alph*h^(1-alph)*gamm/(bet*(1-gamm)-1)));
mu = (c-psii*h^thet*c)^-sig*psii*h^thet/(bet*(1-gamm)-1);
lam = (c-psii*h^thet*c)^-sig*(1+psii*h^thet*gamm/(bet*(1-gamm)-1));

delt1 = alph*(k/h)^(alph-1);
delt2 = delcurv*delt1;

dtarg = dshr - log(log(rbar/rworld)/psii_ln + 1);

eta = lam; 
ETA = lam;

g = 1;

Xss  = [k x i 1 1 1 1 1 1 1 c gdp d];
Yss  = [c x h i 1 lam eta  mu gdp 1 1 1 1 rbar];




%END_EXTRACT_HERE
param.psii = psii;
param.delt1 = delt1;
param.delt2 = delt2;
param.dbar = d;
param.rbar = rbar;
param.dtarg = dtarg;