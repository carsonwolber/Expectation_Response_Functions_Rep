% MODEL_SS - Return the steady state of the model computed analytically
%
% usage:
% 
% [ss, parameters] =model_ss(param)
function [Yss,Xss,param,set] = bll_ss(param,set)

%Upack parameters object
param_unpack


%BEGIN_EXTRACT_HERE


kap = (1-thet*bet)*(1-thet)/thet;
kapw = (1-thetw*bet)*(1-thetw)/(thetw*(1+bet)*(1+zeta*(1+1/muw)));

% r = 1/(bet*gstar^(-1/(1-alph)));  %real rate
% REN = 1/(bet*gstar^(-1/(1-alph))) -1 + delt; %rental rate
% mc = (xi-1)/xi;
% kn = (REN/(mc*alph*gstar))^(1/(alph-1));
% 
% omeg = (1-bet*kapp*gstar^(1/(alph-1)))/(1-kapp*gstar^(1/(alph-1)));  %Coefficient defined in notes
% iy = (1-(1-delt)*gstar^(1/(alph-1)))/(gstar^(alph/(alph-1)))*kn^(1-alph);
% 
% n = ((1-iy-Gshr)/(omeg*mc*(1-alph)))^(-1/(1+1/etta));
% k = kn*n;
% y = gstar^(alph/(alph-1))*k^alph*n^(1-alph);
% c = y*omeg*mc*(1-alph)/(n^(1+1/etta));
% i = iy*y;
% g = Gshr*y;
% w = mc*(1-alph)*gstar^(alph/(alph-1))*k^alph*n^-alph;
% q = 1;
% lam = omeg/c;
% pii = 1;
% ir = r;
% mcbar = mc;
% ikbar = i/k*gstar^(1/(1-alph));
% irbar = ir;
% %Check RC, success!
% %y - c - i - g
% 
% e1 = lam;
% e2 = lam/r;
% e3 = 1;
% e4 = 1;
% 
% v = xi/(xi-1)*mc*y/(1-bet*thet*1^xi);
% pstar = ((1-thet*1^(xi-1))/(1-thet))^(1/(1-xi));
% d = (1-thet)*pstar^(-xi)/(1-thet*1^xi);
 zet =1/(1-psii);
 mu_p = 0.3000;  %where does this come in?
 RkP = 1/bet - (1-delt);
 WAP = ((alph^alph*(1-alph)^(1-alph)/(1+mu_p))*RkP^(-alph))^(1/(1-alph));
 KAN = (WAP/RkP)*(alph/(1-alph));
 IY = delt*(KAN)^(1-alph);
 RkKPY = RkP*KAN^(1-alph);






cybar = (1/zet-IY-RkKPY);
iybar = (IY);
rkpybar = RkKPY;


Xss  = zeros(1,100);
Yss  = zeros(1,100);


bet/((1+bet)*(1-kapw))
%END_EXTRACT_HERE

