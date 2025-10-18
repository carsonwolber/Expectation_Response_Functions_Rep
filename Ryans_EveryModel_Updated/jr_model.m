
%RBC Model code
function [mod,param,set] = jr_model(param,set)

%Name of text files
mod.fname   = 'jr_prog.m';
mod.ss_call = 'jr_ss.m';

%Declare parameters symbols: parameters are values to be estimated, symbols are values that are "fixed"
PARAM = struct2sym(param);
SET   = struct2sym(set);

%Declare Needed Symbols 
syms K   XXL     IL   EPSXN4   EPSXN3   EPSXN2   EPSXN1   GAM  A  GAML CL  GDPL D
syms C   XX   H    I UU LAM   ETA   MU GDP DY DC  DI Q  R

%Declare X and Y vectors
X  = [K   XXL     IL   EPSXN4   EPSXN3   EPSXN2   EPSXN1   GAM A GAML CL  GDPL D];
Y  = [C   XX   H    I UU LAM   ETA   MU  GDP DY DC  DI Q R];

XP = make_prime(X);
YP = make_prime(Y);
make_index([Y,X]);


uc = (C-psii*H^thet*XX)^-sig;


%Model Equations
f(1)     = uc + MU*gamm*(C/XXL)^(gamm-1)*GAML^((gamm-1)/(1-alph)) - LAM;
f(end+1) = uc*psii*thet*H^(thet-1)*XX - LAM*(1-alph)*GAM*A*(UU*K)^alph*H^(-alph);
f(end+1) = uc*psii*H^thet + MU - bet*MU_p*GAM^(-sig/(1-alph))*(C_p/XX)^gamm*(1-gamm)*GAM^(gamm/(1-alph));
f(end+1) = LAM*GAM*alph*A*UU^(alph-1)*K^(alph-1)*H^(1-alph) - ETA*(delt1 + delt2*(UU-1));
f(end+1) = LAM - ETA*(1- phi/2*(I/IL*GAML^(1/(1-alph)) - 1)^2 - phi*(I/IL*GAML^(1/(1-alph)) -1)*I/IL*GAML^(1/(1-alph))) - bet*ETA_p*GAM^(-sig/(1-alph))*phi*(I_p/I*GAM^(1/(1-alph))-1)*(I_p/I*GAM^(1/(1-alph)))^2;
f(end+1) = ETA - bet*GAM^(-sig/(1-alph))*(LAM_p*GAM_p*alph*UU_p^alph*A_p*(K_p/H_p)^(alph-1) + ETA_p*(1-delt0-delt1*(UU_p-1)-delt2/2*(UU_p-1)^2));
f(end+1) = XX - C^gamm*XXL^(1-gamm)*GAML^((gamm-1)/(1-alph));
f(end+1) = GAM*A*(UU*K)^alph*H^(1-alph) + D_p/R - C - I - D;
f(end+1) = K_p*GAM^(1/(1-alph)) - I*(1-phi/2*(I/IL*GAML^(1/(1-alph))-1)^2) - (1-delt0 - delt1*(UU-1)- delt2/2*(UU-1)^2)*K;
f(end+1) = GAM*A*(UU*K)^alph*H^(1-alph) - GDP;
f(end+1) = LAM - bet*(R*LAM_p*GAM^(-sig/(1-alph)));

%OPEN OR CLOSED ECONOMY
if param.open
    f(end+1) = log(R) - log(rworld) - psii_ln*(exp(D/GDP-dtarg) - 1);    
else
    f(end+1) = D_p - dbar;
end

%AUX DEFINITIONS
f(end+1) = DY - GAML^(1/(1-alph))*GDP/GDPL;
f(end+1) = DC - GAML^(1/(1-alph))*C/CL;
f(end+1) = DI - GAML^(1/(1-alph))*I/IL;

f(end+1) = GDPL_p -  GDP;
f(end+1) = CL_p   -  C;
f(end+1) = IL_p   -  I;
f(end+1) = XXL_p  - XX;
f(end+1) = GAML_p - GAM;
f(end+1) = Q - ETA/LAM;

%EXOGENOUS THINGS
f(end+1) = log(EPSXN4_p);
f(end+1) = log(EPSXN3_p) - log(EPSXN4);
f(end+1) = log(EPSXN2_p) - log(EPSXN3);
f(end+1) = log(EPSXN1_p) - log(EPSXN2);
f(end+1) = log(GAM_p)  - rhox*log(GAM) - log(EPSXN1);
f(end+1) = log(A_p)   - .5*log(A) ;

disp(['Neq:  ', num2str(length(f))]);
disp(['Nvar: ', num2str(length(X)+length(Y))]);

%Log-linear approx (Pure linear if log_var = [])
xlog = true(1,length(X));
ylog = true(1,length(Y));
xlog(end) = false;
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
mod.shck(gam_idx-ny  ,1) = siga;
mod.shck(epsxn2_idx-ny,2) = siga;
%Measurement Error (parameters are std-devs in param.m file)
mod.me = [];

%Derivatives using numerical toolbox
mod = anal_deriv(mod);

%Save index variables for use in main_prog
!rm -f rbc_v_idx.mat
neps = size(mod.shck,2);
neq  = nx+ny;


save jr_obj *_idx  mod nx ny neps neq



