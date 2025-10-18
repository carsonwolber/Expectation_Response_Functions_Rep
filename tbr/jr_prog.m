function [f fx fy fxp fyp eta R set dgam_dtheta deta_dtheta dR_dtheta xlag ylag] = model_prog(param, set)

%Assign parameter values to named variables.
sig = param(1);
thet = param(2);
bet = param(3);
rworld = param(4);
hbar = param(5);
gamm = param(6);
alph = param(7);
phi = param(8);
delt0 = param(9);
delcurv = param(10);
delt1 = param(11);
delt2 = param(12);
open = param(13);
dshr = param(14);
psii_ln = param(15);
rhox = param(16);
rhoa = param(17);
rhor = param(18);
sigxn = param(19);
sigxs = param(20);
siga = param(21);
sigr = param(22);
psii = param(23);
kbar = param(24);
dbar = param(25);
abar = param(26);
rbar = param(27);
dtarg = param(28);
e2bar = param(29);
e3bar = param(30);
al = param(31);
gl = param(32);
rl = param(33);
ah = param(34);
gh = param(35);
rh = param(36);

%Assign set values to named variables.
adiff = set(1);
approx_deg = set(2);

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
%Compute Steady State
K= Xss(1);
XXL= Xss(2);
IL= Xss(3);
EPSXN4= Xss(4);
EPSXN3= Xss(5);
EPSXN2= Xss(6);
EPSXN1= Xss(7);
GAM= Xss(8);
A= Xss(9);
GAML= Xss(10);
CL= Xss(11);
GDPL= Xss(12);
D= Xss(13);
C= Yss(1);
XX= Yss(2);
H= Yss(3);
I= Yss(4);
UU= Yss(5);
LAM= Yss(6);
ETA= Yss(7);
MU= Yss(8);
GDP= Yss(9);
DY= Yss(10);
DC= Yss(11);
DI= Yss(12);
Q= Yss(13);
R= Yss(14);
K_p= Xss(1);
XXL_p= Xss(2);
IL_p= Xss(3);
EPSXN4_p= Xss(4);
EPSXN3_p= Xss(5);
EPSXN2_p= Xss(6);
EPSXN1_p= Xss(7);
GAM_p= Xss(8);
A_p= Xss(9);
GAML_p= Xss(10);
CL_p= Xss(11);
GDPL_p= Xss(12);
D_p= Xss(13);
C_p= Yss(1);
XX_p= Yss(2);
H_p= Yss(3);
I_p= Yss(4);
UU_p= Yss(5);
LAM_p= Yss(6);
ETA_p= Yss(7);
MU_p= Yss(8);
GDP_p= Yss(9);
DY_p= Yss(10);
DC_p= Yss(11);
DI_p= Yss(12);
Q_p= Yss(13);
R_p= Yss(14);

%Evaluate F.
f = [1/(C - H^thet*XX*psii)^sig - LAM + (MU*gamm*(C/XXL)^(gamm - 1))/GAML^((gamm - 1)/(alph - 1)), (H^(thet - 1)*XX*psii*thet)/(C - H^thet*XX*psii)^sig + (A*GAM*LAM*(alph - 1)*(K*UU)^alph)/H^alph, MU + (H^thet*psii)/(C - H^thet*XX*psii)^sig + (GAM^(sig/(alph - 1))*MU_p*bet*(C_p/XX)^gamm*(gamm - 1))/GAM^(gamm/(alph - 1)), A*GAM*H^(1 - alph)*K^(alph - 1)*LAM*UU^(alph - 1)*alph - ETA*(delt1 + delt2*(UU - 1)), LAM + ETA*((phi*(I/(GAML^(1/(alph - 1))*IL) - 1)^2)/2 + (I*phi*(I/(GAML^(1/(alph - 1))*IL) - 1))/(GAML^(1/(alph - 1))*IL) - 1) - (ETA_p*GAM^(sig/(alph - 1))*I_p^2*bet*phi*(I_p/(GAM^(1/(alph - 1))*I) - 1))/(GAM^(2/(alph - 1))*I^2), ETA + GAM^(sig/(alph - 1))*bet*(ETA_p*(delt0 + delt1*(UU_p - 1) + (delt2*(UU_p - 1)^2)/2 - 1) - A_p*GAM_p*LAM_p*UU_p^alph*alph*(K_p/H_p)^(alph - 1)), XX - (C^gamm*XXL^(1 - gamm))/GAML^((gamm - 1)/(alph - 1)), D_p/R - D - I - C + A*GAM*H^(1 - alph)*(K*UU)^alph, I*((phi*(I/(GAML^(1/(alph - 1))*IL) - 1)^2)/2 - 1) + K*(delt0 + delt1*(UU - 1) + (delt2*(UU - 1)^2)/2 - 1) + K_p/GAM^(1/(alph - 1)), A*GAM*H^(1 - alph)*(K*UU)^alph - GDP, LAM - GAM^(sig/(alph - 1))*LAM_p*R*bet, D_p - dbar, DY - GDP/(GAML^(1/(alph - 1))*GDPL), DC - C/(CL*GAML^(1/(alph - 1))), DI - I/(GAML^(1/(alph - 1))*IL), GDPL_p - GDP, CL_p - C, IL_p - I, XXL_p - XX, GAML_p - GAM, Q - ETA/LAM, log(EPSXN4_p), log(EPSXN3_p) - log(EPSXN4), log(EPSXN2_p) - log(EPSXN3), log(EPSXN1_p) - log(EPSXN2), log(GAM_p) - log(EPSXN1) - rhox*log(GAM), log(A_p) - log(A)/2];
%Evaluate derivative expressions.
fx = [0, -(C*MU*gamm*(C/XXL)^(gamm - 2)*(gamm - 1))/(GAML^((gamm - 1)/(alph - 1))*XXL), 0, 0, 0, 0, 0, 0, 0, -(GAML*MU*gamm*(C/XXL)^(gamm - 1)*(gamm - 1))/(GAML^((gamm - 1)/(alph - 1) + 1)*(alph - 1)), 0, 0, 0; (A*GAM*K*LAM*UU*alph*(alph - 1)*(K*UU)^(alph - 1))/H^alph, 0, 0, 0, 0, 0, 0, (A*GAM*LAM*(alph - 1)*(K*UU)^alph)/H^alph, (A*GAM*LAM*(alph - 1)*(K*UU)^alph)/H^alph, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, (GAM*GAM^(sig/(alph - 1) - 1)*MU_p*bet*sig*(C_p/XX)^gamm*(gamm - 1))/(GAM^(gamm/(alph - 1))*(alph - 1)) - (GAM*GAM^(sig/(alph - 1))*MU_p*bet*gamm*(C_p/XX)^gamm*(gamm - 1))/(GAM^(gamm/(alph - 1) + 1)*(alph - 1)), 0, 0, 0, 0, 0; A*GAM*H^(1 - alph)*K*K^(alph - 2)*LAM*UU^(alph - 1)*alph*(alph - 1), 0, 0, 0, 0, 0, 0, A*GAM*H^(1 - alph)*K^(alph - 1)*LAM*UU^(alph - 1)*alph, A*GAM*H^(1 - alph)*K^(alph - 1)*LAM*UU^(alph - 1)*alph, 0, 0, 0, 0; 0, 0, -ETA*((I^2*phi)/(GAML^(2/(alph - 1))*IL^2) + (2*I*phi*(I/(GAML^(1/(alph - 1))*IL) - 1))/(GAML^(1/(alph - 1))*IL)), 0, 0, 0, 0, (2*ETA_p*GAM*GAM^(sig/(alph - 1))*I_p^2*bet*phi*(I_p/(GAM^(1/(alph - 1))*I) - 1))/(GAM^(2/(alph - 1) + 1)*I^2*(alph - 1)) + (ETA_p*GAM*GAM^(sig/(alph - 1))*I_p^3*bet*phi)/(GAM^(2/(alph - 1))*GAM^(1/(alph - 1) + 1)*I^3*(alph - 1)) - (ETA_p*GAM*GAM^(sig/(alph - 1) - 1)*I_p^2*bet*phi*sig*(I_p/(GAM^(1/(alph - 1))*I) - 1))/(GAM^(2/(alph - 1))*I^2*(alph - 1)), 0, -ETA*((GAML*I^2*phi)/(GAML^(1/(alph - 1))*GAML^(1/(alph - 1) + 1)*IL^2*(alph - 1)) + (2*GAML*I*phi*(I/(GAML^(1/(alph - 1))*IL) - 1))/(GAML^(1/(alph - 1) + 1)*IL*(alph - 1))), 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, (GAM*GAM^(sig/(alph - 1) - 1)*bet*sig*(ETA_p*(delt0 + delt1*(UU_p - 1) + (delt2*(UU_p - 1)^2)/2 - 1) - A_p*GAM_p*LAM_p*UU_p^alph*alph*(K_p/H_p)^(alph - 1)))/(alph - 1), 0, 0, 0, 0, 0; 0, (C^gamm*XXL*(gamm - 1))/(GAML^((gamm - 1)/(alph - 1))*XXL^gamm), 0, 0, 0, 0, 0, 0, 0, (C^gamm*GAML*XXL^(1 - gamm)*(gamm - 1))/(GAML^((gamm - 1)/(alph - 1) + 1)*(alph - 1)), 0, 0, 0; A*GAM*H^(1 - alph)*K*UU*alph*(K*UU)^(alph - 1), 0, 0, 0, 0, 0, 0, A*GAM*H^(1 - alph)*(K*UU)^alph, A*GAM*H^(1 - alph)*(K*UU)^alph, 0, 0, 0, -1; K*(delt0 + delt1*(UU - 1) + (delt2*(UU - 1)^2)/2 - 1), 0, -(I^2*phi*(I/(GAML^(1/(alph - 1))*IL) - 1))/(GAML^(1/(alph - 1))*IL), 0, 0, 0, 0, -(GAM*K_p)/(GAM^(1/(alph - 1) + 1)*(alph - 1)), 0, -(GAML*I^2*phi*(I/(GAML^(1/(alph - 1))*IL) - 1))/(GAML^(1/(alph - 1) + 1)*IL*(alph - 1)), 0, 0, 0; A*GAM*H^(1 - alph)*K*UU*alph*(K*UU)^(alph - 1), 0, 0, 0, 0, 0, 0, A*GAM*H^(1 - alph)*(K*UU)^alph, A*GAM*H^(1 - alph)*(K*UU)^alph, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -(GAM*GAM^(sig/(alph - 1) - 1)*LAM_p*R*bet*sig)/(alph - 1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, (GAML*GDP)/(GAML^(1/(alph - 1) + 1)*GDPL*(alph - 1)), 0, GDP/(GAML^(1/(alph - 1))*GDPL), 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, (C*GAML)/(CL*GAML^(1/(alph - 1) + 1)*(alph - 1)), C/(CL*GAML^(1/(alph - 1))), 0, 0; 0, 0, I/(GAML^(1/(alph - 1))*IL), 0, 0, 0, 0, 0, 0, (GAML*I)/(GAML^(1/(alph - 1) + 1)*IL*(alph - 1)), 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -GAM, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -1, -rhox, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -1/2, 0, 0, 0, 0];
fy = [(C*MU*gamm*(C/XXL)^(gamm - 2)*(gamm - 1))/(GAML^((gamm - 1)/(alph - 1))*XXL) - (C*sig)/(C - H^thet*XX*psii)^(sig + 1), (H^thet*XX*psii*sig)/(C - H^thet*XX*psii)^(sig + 1), (H*H^(thet - 1)*XX*psii*sig*thet)/(C - H^thet*XX*psii)^(sig + 1), 0, 0, -LAM, 0, (MU*gamm*(C/XXL)^(gamm - 1))/GAML^((gamm - 1)/(alph - 1)), 0, 0, 0, 0, 0, 0; -(C*H^(thet - 1)*XX*psii*sig*thet)/(C - H^thet*XX*psii)^(sig + 1), (H^(thet - 1)*XX*psii*thet)/(C - H^thet*XX*psii)^sig + (H^thet*H^(thet - 1)*XX^2*psii^2*sig*thet)/(C - H^thet*XX*psii)^(sig + 1), (H*H^(thet - 2)*XX*psii*thet*(thet - 1))/(C - H^thet*XX*psii)^sig + (H*H^(2*thet - 2)*XX^2*psii^2*sig*thet^2)/(C - H^thet*XX*psii)^(sig + 1) - (A*GAM*H*LAM*alph*(alph - 1)*(K*UU)^alph)/H^(alph + 1), 0, (A*GAM*K*LAM*UU*alph*(alph - 1)*(K*UU)^(alph - 1))/H^alph, (A*GAM*LAM*(alph - 1)*(K*UU)^alph)/H^alph, 0, 0, 0, 0, 0, 0, 0, 0; -(C*H^thet*psii*sig)/(C - H^thet*XX*psii)^(sig + 1), (H^(2*thet)*XX*psii^2*sig)/(C - H^thet*XX*psii)^(sig + 1) - (C_p*GAM^(sig/(alph - 1))*MU_p*bet*gamm*(C_p/XX)^(gamm - 1)*(gamm - 1))/(GAM^(gamm/(alph - 1))*XX), (H*H^(thet - 1)*psii*thet)/(C - H^thet*XX*psii)^sig + (H*H^thet*H^(thet - 1)*XX*psii^2*sig*thet)/(C - H^thet*XX*psii)^(sig + 1), 0, 0, 0, 0, MU, 0, 0, 0, 0, 0, 0; 0, 0, -(A*GAM*H*K^(alph - 1)*LAM*UU^(alph - 1)*alph*(alph - 1))/H^alph, 0, A*GAM*H^(1 - alph)*K^(alph - 1)*LAM*UU*UU^(alph - 2)*alph*(alph - 1) - ETA*UU*delt2, A*GAM*H^(1 - alph)*K^(alph - 1)*LAM*UU^(alph - 1)*alph, -ETA*(delt1 + delt2*(UU - 1)), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, ETA*((I^2*phi)/(GAML^(2/(alph - 1))*IL^2) + (2*I*phi*(I/(GAML^(1/(alph - 1))*IL) - 1))/(GAML^(1/(alph - 1))*IL)) + (ETA_p*GAM^(sig/(alph - 1))*I_p^3*bet*phi)/(GAM^(3/(alph - 1))*I^3) + (2*ETA_p*GAM^(sig/(alph - 1))*I_p^2*bet*phi*(I_p/(GAM^(1/(alph - 1))*I) - 1))/(GAM^(2/(alph - 1))*I^2), 0, LAM, ETA*((phi*(I/(GAML^(1/(alph - 1))*IL) - 1)^2)/2 + (I*phi*(I/(GAML^(1/(alph - 1))*IL) - 1))/(GAML^(1/(alph - 1))*IL) - 1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, ETA, 0, 0, 0, 0, 0, 0, 0; -(C*C^(gamm - 1)*XXL^(1 - gamm)*gamm)/GAML^((gamm - 1)/(alph - 1)), XX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -C, 0, -(A*GAM*H*(alph - 1)*(K*UU)^alph)/H^alph, -I, A*GAM*H^(1 - alph)*K*UU*alph*(K*UU)^(alph - 1), 0, 0, 0, 0, 0, 0, 0, 0, -D_p/R; 0, 0, 0, I*((phi*(I/(GAML^(1/(alph - 1))*IL) - 1)^2)/2 - 1) + (I^2*phi*(I/(GAML^(1/(alph - 1))*IL) - 1))/(GAML^(1/(alph - 1))*IL), K*(UU*delt1 + UU*delt2*(UU - 1)), 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -(A*GAM*H*(alph - 1)*(K*UU)^alph)/H^alph, 0, A*GAM*H^(1 - alph)*K*UU*alph*(K*UU)^(alph - 1), 0, 0, 0, -GDP, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, LAM, 0, 0, 0, 0, 0, 0, 0, -GAM^(sig/(alph - 1))*LAM_p*R*bet; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -GDP/(GAML^(1/(alph - 1))*GDPL), DY, 0, 0, 0, 0; -C/(CL*GAML^(1/(alph - 1))), 0, 0, 0, 0, 0, 0, 0, 0, 0, DC, 0, 0, 0; 0, 0, 0, -I/(GAML^(1/(alph - 1))*IL), 0, 0, 0, 0, 0, 0, 0, DI, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -GDP, 0, 0, 0, 0, 0; -C, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -I, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -XX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, ETA/LAM, -ETA/LAM, 0, 0, 0, 0, 0, Q, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
fxp = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -(A_p*GAM^(sig/(alph - 1))*GAM_p*K_p*LAM_p*UU_p^alph*alph*bet*(K_p/H_p)^(alph - 2)*(alph - 1))/H_p, 0, 0, 0, 0, 0, 0, -A_p*GAM^(sig/(alph - 1))*GAM_p*LAM_p*UU_p^alph*alph*bet*(K_p/H_p)^(alph - 1), -A_p*GAM^(sig/(alph - 1))*GAM_p*LAM_p*UU_p^alph*alph*bet*(K_p/H_p)^(alph - 1), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/R; K_p/GAM^(1/(alph - 1)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, GDPL_p, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, CL_p, 0, 0; 0, 0, IL_p, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, XXL_p, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, GAML_p, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];
fyp = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; (C_p*GAM^(sig/(alph - 1))*MU_p*bet*gamm*(C_p/XX)^(gamm - 1)*(gamm - 1))/(GAM^(gamm/(alph - 1))*XX), 0, 0, 0, 0, 0, 0, (GAM^(sig/(alph - 1))*MU_p*bet*(C_p/XX)^gamm*(gamm - 1))/GAM^(gamm/(alph - 1)), 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, - (ETA_p*GAM^(sig/(alph - 1))*I_p^3*bet*phi)/(GAM^(3/(alph - 1))*I^3) - (2*ETA_p*GAM^(sig/(alph - 1))*I_p^2*bet*phi*(I_p/(GAM^(1/(alph - 1))*I) - 1))/(GAM^(2/(alph - 1))*I^2), 0, 0, -(ETA_p*GAM^(sig/(alph - 1))*I_p^2*bet*phi*(I_p/(GAM^(1/(alph - 1))*I) - 1))/(GAM^(2/(alph - 1))*I^2), 0, 0, 0, 0, 0, 0, 0; 0, 0, (A_p*GAM^(sig/(alph - 1))*GAM_p*K_p*LAM_p*UU_p^alph*alph*bet*(K_p/H_p)^(alph - 2)*(alph - 1))/H_p, 0, GAM^(sig/(alph - 1))*bet*(ETA_p*(UU_p*delt1 + UU_p*delt2*(UU_p - 1)) - A_p*GAM_p*LAM_p*UU_p*UU_p^(alph - 1)*alph^2*(K_p/H_p)^(alph - 1)), -A_p*GAM^(sig/(alph - 1))*GAM_p*LAM_p*UU_p^alph*alph*bet*(K_p/H_p)^(alph - 1), ETA_p*GAM^(sig/(alph - 1))*bet*(delt0 + delt1*(UU_p - 1) + (delt2*(UU_p - 1)^2)/2 - 1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -GAM^(sig/(alph - 1))*LAM_p*R*bet, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

eta = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, siga; 0, 0; siga, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0];
R = [];
