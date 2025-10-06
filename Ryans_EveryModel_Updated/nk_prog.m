function [f fx fy fxp fyp eta R set dgam_dtheta deta_dtheta dR_dtheta xlag ylag] = model_prog(param, set)

%Assign parameter values to named variables.
rhoga = param(1);
sigga = param(2);
sigea = param(3);
siges = param(4);

%Assign set values to named variables.
adiff = set(1);
approx_deg = set(2);
bet = set(3);
alph = set(4);
thet = set(5);
etap = set(6);
zeta = set(7);
rhoe = set(8);
phipi = set(9);
phiy = set(10);
phidy = set(11);
rhoi = set(12);
gstar = set(13);

%BEGIN_EXTRACT_HERE




%Y and X vectors with SS values
Yss = zeros(1,100);
Xss = zeros(1,100);

ss = [Yss Xss];

%END_EXTRACT_HERE
%Compute Steady State
IRL= Xss(1);
GAPL= Xss(2);
EI= Xss(3);
AL= Xss(4);
GAM= Xss(5);
A= Xss(6);
GAP= Yss(1);
IR= Yss(2);
PI= Yss(3);
DY= Yss(4);
IRL_p= Xss(1);
GAPL_p= Xss(2);
EI_p= Xss(3);
AL_p= Xss(4);
GAM_p= Xss(5);
A_p= Xss(6);
GAP_p= Yss(1);
IR_p= Yss(2);
PI_p= Yss(3);
DY_p= Yss(4);

%Evaluate F.
f = [PI - PI_p*bet - (GAP*(bet*thet - 1)*(thet - 1)*(zeta + 1))/(thet*(alph*(etap - 1) + 1)), GAP - GAP_p + IR - PI_p + (alph - 1)*(A_p - A + GAM_p), IR - EI - IRL*rhoi + (rhoi - 1)*(DY*phidy + GAP*phiy + PI*phipi), IRL_p - IR, GAPL_p - GAP, DY - GAP + GAPL + (alph - 1)*(A - AL + GAM), AL_p - A, EI_p, GAM_p, A_p];
%Evaluate derivative expressions.
fx = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1 - alph; -rhoi, 0, -1, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 1 - alph, alph - 1, alph - 1; 0, 0, 0, 0, 0, -1; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0];
fy = [-((bet*thet - 1)*(thet - 1)*(zeta + 1))/(thet*(alph*(etap - 1) + 1)), 0, 1, 0; 1, 1, 0, 0; phiy*(rhoi - 1), 1, phipi*(rhoi - 1), phidy*(rhoi - 1); 0, -1, 0, 0; -1, 0, 0, 0; -1, 0, 0, 1; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0];
fxp = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, alph - 1, alph - 1; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0; 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, 1];
fyp = [0, 0, -bet, 0; -1, 0, -1, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0];

eta = [1; 0; 0; 0; 0; 0];
R = [];
