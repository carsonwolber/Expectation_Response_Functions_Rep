%RBC Model code
 

function [mod,param,set] = nk_model(param,set)


%Name of text files
mod.fname   = 'nk_prog.m';
mod.ss_call = 'nk_ss.m';

%Declare parameters symbols: parameters are values to be estimated, symbols are values that are "fixed"
PARAM = struct2sym(param);
SET   = struct2sym(set);

%Declare Needed Symbols 
syms IRL GAPL EI GAM
syms GAP PI DY A AL
IR = sym('IR');
%Declare X and Y vectors
X  = [IRL GAPL EI AL GAM A];
Y  = [GAP IR PI DY];
XP = make_prime(X);
YP = make_prime(Y);
make_index([Y,X]);


kap = (1-bet*thet)*(1-thet)/thet*(1+zeta)/(1+alph*(etap-1));

%Model Equations
f(1)     = PI - kap*GAP - bet*PI_p;
f(end+1) = GAP + IR - (GAP_p + PI_p + (1-alph)*(GAM_p + A_p-A));
f(end+1) = IR - rhoi*IRL - (1-rhoi)*(phipi*PI + phiy*GAP + phidy*DY)  - EI;


%AUX DEFINITIONS
f(end+1) = IRL_p - IR;
f(end+1) = GAPL_p  - GAP;
f(end+1) = DY    - (GAP - GAPL + (1-alph)*(GAM + A-AL));
f(end+1) = AL_p - A;
f(end+1) = EI_p;
f(end+1) = GAM_p;
f(end+1) = A_p;



disp(['Neq:  ', num2str(length(f))]);
disp(['Nvar: ', num2str(length(X)+length(Y))]);

%Log-linear approx (Pure linear if log_var = [])
xlog = false(1,length(X));
ylog = false(1,length(Y));
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
mod.shck(1) = 1;

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


save nk_obj *_idx  mod nx ny neps neq



