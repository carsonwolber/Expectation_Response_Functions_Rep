% MODEL_SS - Return the steady state of the model computed analytically
%
% usage:
% 
% [ss, parameters] =model_ss(param)


function [ss,param,set] = rbc_ss(param,set)

%Upack parameters object
param_unpack


%BEGIN_EXTRACT_HERE


%Steady state, use closed form expressions for the ss values.
r = (1/(bet*gam^(-1/(1-alph))))-1+del;
kh = (r/(alph*gam))^(1/(alph-1));
w = (1-alph)*(r/alph)^(alph/(alph-1));
ik = (1-(1-del)/(gam^(1/(1-alph))));
c = w/chi;
k = c/(gam^(alph/(alph-1))*kh^(alph-1) - ik);
i = ik*k;
h = kh^-1*k;
gdp = c+i;

%Y and X vectors with SS values
Yss = [gdp c h w r i gam^(1/(1-alph))*ones(1,3)];
Xss = [gam 1 k c i gdp];


ss = [Yss Xss];

%END_EXTRACT_HERE
