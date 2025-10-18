%RBC Model code
function [mod,param,set] = bll_model(param,set)

%Name of text files
mod.fname   = 'bll_prog.m';
mod.ss_call = 'bll_ss.m';

%Declare parameters symbols: parameters are values to be estimated, symbols are values that are "fixed"
PARAM = struct2sym(param);
SET   = struct2sym(set);

%Declare Needed Symbols
syms KBAR CL IL IRL WL D GAM MPT MWT G GDPL DY DC DI GAM Q EPSP EPSW SIGNAL RREAL
syms C I GDP N K U LAM PHI RK W PI GAP IRT DY DW DN NL DX Z ZL V DXT ZT ZLT VT
syms DX DZ ETTA V VL DXT DZT ETTAT VT VLT C_p LAMPT PIPT  PIPT PHIPT RKPT IPT WPT A AL

syms DELAT MP3T MP2T MP1T MMT ML1T ML2T XIT XILT EXT EXLT EXL2T DELA MP3 MP2 MP1 MM ML1 ML2 XI XIL EX EXL EXL2 DELAT MP3T MP2T MP1T MT ML1T ML2T XIT XILT EXT EXLT EXL2T
%Declare X and Y vectors
X  = [KBAR CL IL IRL WL NL GDPL AL D MPT MWT G Q EPSP EPSW  GAM A];
Y  = [DY DC DI DW DN PI C I GDP N K U LAM PHI RK W  GAP IRT RREAL ...
    LAMPT PIPT PHIPT RKPT IPT WPT]; %Variables that appear as expectations
XP = make_prime(X);
YP = make_prime(Y);
make_index([Y,X]);

%Model Equations
GAMM   = GAM   + A   - AL;
GAMM_p = GAM_p + A_p - A;

f(1)     = (gam-h*bet)*(gam-h)*LAM - h*bet*gam*C_p + (gam^2+h^2*bet)*C - h*gam*CL - h*bet*gam*GAMM_p + h*gam*GAMM;
f(end+1) = LAM-IRT-LAMPT+GAMM_p+PIPT;
f(end+1) = PHI-(1-delt)*bet*gam^-1*(PHIPT-GAMM_p)-(1-(1-delt)*bet*gam^-1)*(LAMPT - GAMM_p + RKPT);
f(end+1) = LAM - PHI - D + chi*gam^2*(I-IL+GAMM) - bet*chi*gam^2*(IPT - I + GAMM_p);
f(end+1) = RK - xi*U;
f(end+1) = GAP - alph*RK - (1-alph)*W;
f(end+1) = RK - W + K - N;
f(end+1) = K - U - KBAR + GAMM;
f(end+1) = KBAR_p - (1-delt)*gam^-1*(KBAR-GAMM) - (1-(1-delt)*gam^-1)*(D+I);
f(end+1) = GDP - alph*K - (1-alph)*N;
f(end+1) = (1-psii)*GDP - cybar*C - iybar*I - rkpybar*U - (1-psii)*G; 
f(end+1) = PI - bet*PIPT - kap*GAP - MPT;
f(end+1) = W - 1/(1+bet)*WL - bet/(1+bet)*(WPT) + 1/(1+bet)*(PI+GAMM) - bet/(1+bet)*(PIPT + GAMM_p) + kapw*(W - zeta*N + LAM) - kapw*MWT;
f(end+1) = IRT - rhoi*IRL - (1-rhoi)*(phipi*PI + phiy*GDP + phidy*DY) - Q;


%Impose alternative PLM
f(end+1) = LAMPT - LAM_p;
f(end+1) = PIPT  - PI_p;

f(end+1) = PHIPT - PHI_p  ;
f(end+1) = RKPT  - RK_p   ;
f(end+1) = IPT   - I_p    ;
f(end+1) = WPT   - W_p    ;

%Real interest rate
f(end+1) = RREAL - (IRT - PIPT);

%AUX DEFINITIONS
f(end+1) = CL_p  - C;
f(end+1) = IRL_p - IRT;
f(end+1) = IL_p - I;
f(end+1) = WL_p  - W;
f(end+1) = NL_p  - N;
f(end+1) = GDPL_p - GDP;
f(end+1) = DY - GDP + GDPL - GAMM;
f(end+1) = DC - C + CL - GAMM;
f(end+1) = DI - I + IL - GAMM;
f(end+1) = DW - W + WL - GAMM;
f(end+1) = DN - N + NL;
f(end+1) = AL_p - A;

%EXOGN PROCESS
f(end+1) = EPSW_p;
f(end+1) = EPSP_p;


f(end+1) = D_p - rhod*D;
f(end+1) = MWT_p - rhow*MWT - phiw*EPSW;
f(end+1) = MPT_p - rhop*MPT - phip*EPSP;
f(end+1) = Q_p-rhoq*Q;
f(end+1) = G_p-rhog*G;
f(end+1) = GAM_p;
f(end+1) = A_p;

disp(['Neq:  ', num2str(length(f))]);
disp(['Nvar: ', num2str(length(X)+length(Y))]);

%Log-linear approx (Pure linear if log_var = [])
xlog = false(1,length(X));
ylog = false(1,length(Y));
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

%Measurement Error (parameters are std-devs in param.m file)
mod.me = [];

%Derivatives using numerical toolbox
mod = anal_deriv(mod);

%Save index variables for use in main_prog
neps = size(mod.shck,2);
neq  = nx+ny;
save bll_obj *_idx  mod nx ny neps neq



