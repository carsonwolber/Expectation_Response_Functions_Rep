%RBC Model code
 

function [mod,param,set] = rbc_model(param,set)

%Name of text files
mod.fname   = 'rbc_prog.m';
mod.ss_call = 'rbc_ss.m';

%Declare parameters symbols: parameters are values to be estimated, symbols are values that are "fixed"
PARAM = struct2sym(param);
SET   = struct2sym(set);

%Declare Needed Symbols
syms GAM K A 
syms GDP C H W R  I  CL IL DC DI YL DY

%Declare X and Y vectors
X  = [GAM A  K CL IL YL];
Y  = [GDP C   H   W   R   I DC DI DY]; 
XP = make_prime(X);
YP = make_prime(Y);
make_index([Y,X]);


%Model Equations
f(1)     = C-W/chi;
f(end+1) = 1-bet*C/C_p*GAM_p^(-1/(1-alph))*(R_p+1-del);
f(end+1) = R - alph*GAM*A*(K/H)^(alph-1);
f(end+1) = W - (1-alph)*GAM^(alph/(alph-1))*A*(K/H)^alph;
f(end+1) = (1-del)*K + I*GAM^(1/(1-alph)) - K_p*GAM^(1/(1-alph));
f(end+1) = C+I- GAM^(alph/(alph-1))*A*K^alph*H^(1-alph);

f(end+1) = GDP - GAM^(alph/(alph-1))*A*K^alph*H^(1-alph);
f(end+1) = CL_p - C;
f(end+1) = IL_p - I;
f(end+1) = YL_p - GDP;
f(end+1) = DC - C/CL  *GAM^(1/(1-alph));
f(end+1) = DI - I/IL  *GAM^(1/(1-alph));
f(end+1) = DY - GDP/YL*GAM^(1/(1-alph));
f(end+1) = log(GAM_p/gam) - rho*log(GAM/gam);
f(end+1) = log(A_p) - .9*log(A);

disp(['Neq:  ', num2str(length(f))]);
disp(['Nvar: ', num2str(length(X)+length(Y))]);

%Log-linear approx (Pure linear if log_var = [])
xlog = true(1,length(X));
ylog = true(1,length(Y));
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
mod.shck = sym(zeros(nx,1));
mod.shck(1) = siga;

%Measurement Error (parameters are std-devs in param.m file)
mod.me = [];

%Derivatives using numerical toolbox
mod = anal_deriv(mod);

%Save index variables for use in main_prog
!rm -f rbc_v_idx.mat
nx   = length(X);
ny   = length(Y);
neps = size(mod.shck,2);
neq  = nx+ny;


save rbc_obj *_idx  mod nx ny neps neq



