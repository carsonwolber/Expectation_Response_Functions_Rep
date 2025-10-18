%RBC Model code
function [mod,param,set] = bs_model(param,set)

%Name of text files
mod.fname   = 'bs_prog.m';
mod.ss_call = 'bs_ss.m';

%Declare parameters symbols: parameters are values to be estimated, symbols are values that are "fixed"
PARAM = struct2sym(param);
SET   = struct2sym(set);

%Declare Needed Symbols 
syms K DA GG GSHR CL IRL YL EI ES MU1 GGL MU1T GGT  XX XXL XXLL V EV XXT XXLT XXLLT VT EVT DAT
syms C I GDP GOV H W REN R PI LAM Q MC DY E1 E2 E3 E4 IL DC DI SIGNAL DN HL DAPT
IR = sym('IR');
%Declare X and Y vectors
X  = [K GSHR CL IL IRL YL HL EI ES DA]
Y  = [C I GDP GOV H W REN R IR PI LAM Q MC DY DC DI DN E1 E2 E3 E4 ];
XP = make_prime(X);
YP = make_prime(Y);
make_index([Y,X]);

%Adjustment cost
PHI = @(x) x - gam/2*(x-ikbar)^2;
PHIP = @(x) 1 - gam*(x-ikbar);

%Model Equations
f(1)     = LAM - E1;
f(end+1) = LAM/R - E2;
f(end+1) = H^(1/etta) - LAM*W;
f(end+1) = W - MC*(1-alph)*DA^(alph/(alph-1))*K^alph*H^(-alph);
f(end+1) = REN - MC*( alph )*DA                *K^(alph-1)*H^(1-alph);
f(end+1) = PI - (mcbar^((1-thet)/thet*(1-bet*thet)))^-1*MC^((1-thet)/thet*(1-bet*thet))*E3^bet;
f(end+1) = Q*PHIP(I/K*DA^(1/(1-alph)))-1;
f(end+1) = Q - E4;
f(end+1) = GDP - DA^(alph/(alph-1))*K^alph*H^(1-alph);
f(end+1) = GDP - C - I - GOV;
f(end+1) = R - IR/E3;
f(end+1) = log(IR/irbar) - rhoi*log(IRL/irbar) - (1-rhoi)*(phipi*log(PI) + phiy*log(GDP/gdp) + phidy*log(DY/gstar^(1/(1-alph)))) - log(EI);
f(end+1) = K_p - (1-delt)*DA^(1/(alph-1))*K - I;

%Define expectations terms
f(end+1) = E1 - (C-kapp*CL*DA^(1/(alph-1)))^-1 + bet*kapp*(C_p*DA_p^(1/(1-alph)) - kapp*C)^-1;
f(end+1) = E2 - bet*LAM_p*DA_p^(1/(alph-1));
f(end+1) = E3 - PI_p;
f(end+1) = E4 - bet*(LAM_p/LAM*DA_p^(1/(alph-1))*(REN_p + Q_p*((1-delt) + PHI(I_p/K_p*DA_p^(1/(1-alph))) - I_p/K_p*DA_p^(1/(1-alph))*PHIP(I_p/K_p*DA_p^(1/(1-alph))))));




%AUX DEFINITIONS
f(end+1) = CL_p  - C;
f(end+1) = IRL_p - IR;
f(end+1) = YL_p  - GDP;
f(end+1) = IL_p  - I;
f(end+1) = HL_p  - H;
f(end+1) = DY    - GDP/YL*DA^(1/(1-alph));
f(end+1) = DC    - C/CL*DA^(1/(1-alph));
f(end+1) = DI    - I/IL*DA^(1/(1-alph));
f(end+1) = DN    - H/HL;
f(end+1) = GOV/GDP - GSHR;
f(end+1) = log(GSHR_p/Gshr) - rhog*log(GSHR/Gshr);

%EXOG PROCESSES
f(end+1) = log(EI_p);
f(end+1) = log(ES_p);
f(end+1) = log(DA_p/gstar);

disp(['Neq:  ', num2str(length(f))]);
disp(['Nvar: ', num2str(length(X)+length(Y))]);

%Log-linear approx (Pure linear if log_var = [])
xlog = true(1,length(X));
ylog = true(1,length(Y));
%xlog(end) = false;
log_var = [X(xlog) Y(ylog) XP(xlog) YP(ylog)];


mod.f  = subs(f, log_var, exp(log_var));
mod.X  = X;
mod.XP = XP;
mod.Y  = Y;
mod.YP = YP;
mod.PARAM = PARAM;
mod.param = param;
mod.SET = SET;
mod.set = set;
mod.adiff = set.adiff; %Include anaylytical derivatives?
mod.xlog = xlog;
mod.ylog = ylog;


%Standard Errors
nx = length(X);
ny = length(Y);
mod.shck = sym(zeros(nx,2));
mod.shck(da_idx-ny  ,1)  = sigga;
mod.shck(ei_idx-ny  ,2)  = 1;

%Measurement Error (parameters are std-devs in param.m file)
mod.me = [];

%Derivatives using numerical toolbox
mod = anal_deriv(mod);

%Save index variables for use in main_prog
!rm -f rbc_v_idx.mat
neps = size(mod.shck,2);
neq  = nx+ny;
save bs_obj *_idx  mod nx ny neps neq



