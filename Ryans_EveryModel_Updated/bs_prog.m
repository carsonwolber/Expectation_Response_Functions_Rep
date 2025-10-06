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
delt = set(5);
gstar = set(6);
Gshr = set(7);
rhog = set(8);
sigg = set(9);
rhoe = set(10);
phipi = set(11);
phiy = set(12);
phidy = set(13);
rhoi = set(14);
kapp = set(15);
gam = set(16);
etta = set(17);
sigi = set(18);
xi = set(19);
thet = set(20);
ybar = set(21);
mcbar = set(22);
ikbar = set(23);
irbar = set(24);
gdp = set(25);

%BEGIN_EXTRACT_HERE
r = 1/(bet*gstar^(-1/(1-alph)))  %real rate
REN = 1/(bet*gstar^(-1/(1-alph))) -1 + delt; %rental rate
mc = (xi-1)/xi;
kn = (REN/(mc*alph*gstar))^(1/(alph-1));

omeg = (1-bet*kapp*gstar^(1/(alph-1)))/(1-kapp*gstar^(1/(alph-1)));  %Coefficient defined in notes
iy = (1-(1-delt)*gstar^(1/(alph-1)))/(gstar^(alph/(alph-1)))*kn^(1-alph);

n = ((1-iy-Gshr)/(omeg*mc*(1-alph)))^(-1/(1+1/etta));
k = kn*n;
y = gstar^(alph/(alph-1))*k^alph*n^(1-alph);
c = y*omeg*mc*(1-alph)/(n^(1+1/etta));
i = iy*y;
g = Gshr*y;
w = mc*(1-alph)*gstar^(alph/(alph-1))*k^alph*n^-alph;
q = 1;
lam = omeg/c;
pii = 1;
ir = r;
mcbar = mc;
ikbar = i/k*gstar^(1/(1-alph));
irbar = ir;
%Check RC, success!
%y - c - i - g
gdp = y;

e1 = lam;
e2 = lam/r;
e3 = 1;
e4 = 1;

v = xi/(xi-1)*mc*y/(1-bet*thet*1^xi);
pstar = ((1-thet*1^(xi-1))/(1-thet))^(1/(1-xi));
d = (1-thet)*pstar^(-xi)/(1-thet*1^xi);

Xss  = [k  Gshr c i ir y n 1 1 gstar];
Yss  = [c i y g n w REN r ir pii lam q mc gstar^(1/(1-alph)) gstar^(1/(1-alph)) gstar^(1/(1-alph)) 1 e1 e2 e3 e4];

ss = [Yss Xss];


%END_EXTRACT_HERE
%Compute Steady State
K= Xss(1);
GSHR= Xss(2);
CL= Xss(3);
IL= Xss(4);
IRL= Xss(5);
YL= Xss(6);
HL= Xss(7);
EI= Xss(8);
ES= Xss(9);
DA= Xss(10);
C= Yss(1);
I= Yss(2);
GDP= Yss(3);
GOV= Yss(4);
H= Yss(5);
W= Yss(6);
REN= Yss(7);
R= Yss(8);
IR= Yss(9);
PI= Yss(10);
LAM= Yss(11);
Q= Yss(12);
MC= Yss(13);
DY= Yss(14);
DC= Yss(15);
DI= Yss(16);
DN= Yss(17);
E1= Yss(18);
E2= Yss(19);
E3= Yss(20);
E4= Yss(21);
K_p= Xss(1);
GSHR_p= Xss(2);
CL_p= Xss(3);
IL_p= Xss(4);
IRL_p= Xss(5);
YL_p= Xss(6);
HL_p= Xss(7);
EI_p= Xss(8);
ES_p= Xss(9);
DA_p= Xss(10);
C_p= Yss(1);
I_p= Yss(2);
GDP_p= Yss(3);
GOV_p= Yss(4);
H_p= Yss(5);
W_p= Yss(6);
REN_p= Yss(7);
R_p= Yss(8);
IR_p= Yss(9);
PI_p= Yss(10);
LAM_p= Yss(11);
Q_p= Yss(12);
MC_p= Yss(13);
DY_p= Yss(14);
DC_p= Yss(15);
DI_p= Yss(16);
DN_p= Yss(17);
E1_p= Yss(18);
E2_p= Yss(19);
E3_p= Yss(20);
E4_p= Yss(21);

%Evaluate F.
f = [LAM - E1, LAM/R - E2, H^(1/etta) - LAM*W, W + (DA^(alph/(alph - 1))*K^alph*MC*(alph - 1))/H^alph, REN - DA*H^(1 - alph)*K^(alph - 1)*MC*alph, PI - (E3^bet*MC^(((bet*thet - 1)*(thet - 1))/thet))/mcbar^(((bet*thet - 1)*(thet - 1))/thet), Q*(gam*(ikbar - I/(DA^(1/(alph - 1))*K)) + 1) - 1, Q - E4, GDP - DA^(alph/(alph - 1))*H^(1 - alph)*K^alph, GDP - C - GOV - I, R - IR/E3, log(IR/irbar) - log(EI) + (rhoi - 1)*(phipi*log(PI) + phiy*log(GDP/gdp) + phidy*log(DY*gstar^(1/(alph - 1)))) - rhoi*log(IRL/irbar), K_p - I + DA^(1/(alph - 1))*K*(delt - 1), E1 - 1/(C - CL*DA^(1/(alph - 1))*kapp) - (bet*kapp)/(C*kapp - C_p/DA_p^(1/(alph - 1))), E2 - DA_p^(1/(alph - 1))*LAM_p*bet, E3 - PI_p, E4 - (DA_p^(1/(alph - 1))*LAM_p*bet*(REN_p - Q_p*(delt + (gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p))^2)/2 - I_p/(DA_p^(1/(alph - 1))*K_p) + (I_p*(gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)) + 1))/(DA_p^(1/(alph - 1))*K_p) - 1)))/LAM, CL_p - C, IRL_p - IR, YL_p - GDP, IL_p - I, HL_p - H, DY - GDP/(DA^(1/(alph - 1))*YL), DC - C/(CL*DA^(1/(alph - 1))), DI - I/(DA^(1/(alph - 1))*IL), DN - H/HL, GOV/GDP - GSHR, log(GSHR_p/Gshr) - rhog*log(GSHR/Gshr), log(EI_p), log(ES_p), log(DA_p/gstar)];
%Evaluate derivative expressions.
fx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; (DA^(alph/(alph - 1))*K*K^(alph - 1)*MC*alph*(alph - 1))/H^alph, 0, 0, 0, 0, 0, 0, 0, 0, (DA*DA^(alph/(alph - 1) - 1)*K^alph*MC*alph)/H^alph; -DA*H^(1 - alph)*K*K^(alph - 2)*MC*alph*(alph - 1), 0, 0, 0, 0, 0, 0, 0, 0, -DA*H^(1 - alph)*K^(alph - 1)*MC*alph; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; (I*Q*gam)/(DA^(1/(alph - 1))*K), 0, 0, 0, 0, 0, 0, 0, 0, (DA*I*Q*gam)/(DA^(1/(alph - 1) + 1)*K*(alph - 1)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -DA^(alph/(alph - 1))*H^(1 - alph)*K*K^(alph - 1)*alph, 0, 0, 0, 0, 0, 0, 0, 0, -(DA*DA^(alph/(alph - 1) - 1)*H^(1 - alph)*K^alph*alph)/(alph - 1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -rhoi, 0, 0, -1, 0, 0; DA^(1/(alph - 1))*K*(delt - 1), 0, 0, 0, 0, 0, 0, 0, 0, (DA*DA^(1/(alph - 1) - 1)*K*(delt - 1))/(alph - 1); 0, 0, -(CL*DA^(1/(alph - 1))*kapp)/(C - CL*DA^(1/(alph - 1))*kapp)^2, 0, 0, 0, 0, 0, 0, -(CL*DA*DA^(1/(alph - 1) - 1)*kapp)/((C - CL*DA^(1/(alph - 1))*kapp)^2*(alph - 1)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, GDP/(DA^(1/(alph - 1))*YL), 0, 0, 0, (DA*GDP)/(DA^(1/(alph - 1) + 1)*YL*(alph - 1)); 0, 0, C/(CL*DA^(1/(alph - 1))), 0, 0, 0, 0, 0, 0, (C*DA)/(CL*DA^(1/(alph - 1) + 1)*(alph - 1)); 0, 0, 0, I/(DA^(1/(alph - 1))*IL), 0, 0, 0, 0, 0, (DA*I)/(DA^(1/(alph - 1) + 1)*IL*(alph - 1)); 0, 0, 0, 0, 0, 0, H/HL, 0, 0, 0; 0, -GSHR, 0, 0, 0, 0, 0, 0, 0, 0; 0, -rhog, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
fy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, LAM, 0, 0, 0, 0, 0, 0, -E1, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -LAM/R, 0, 0, LAM/R, 0, 0, 0, 0, 0, 0, 0, -E2, 0, 0; 0, 0, 0, 0, (H*H^(1/etta - 1))/etta, -LAM*W, 0, 0, 0, 0, -LAM*W, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -(DA^(alph/(alph - 1))*H*K^alph*MC*alph*(alph - 1))/H^(alph + 1), W, 0, 0, 0, 0, 0, 0, (DA^(alph/(alph - 1))*K^alph*MC*(alph - 1))/H^alph, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, (DA*H*K^(alph - 1)*MC*alph*(alph - 1))/H^alph, 0, REN, 0, 0, 0, 0, 0, -DA*H^(1 - alph)*K^(alph - 1)*MC*alph, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, PI, 0, 0, -(E3^bet*MC*MC^(((bet*thet - 1)*(thet - 1))/thet - 1)*(bet*thet - 1)*(thet - 1))/(mcbar^(((bet*thet - 1)*(thet - 1))/thet)*thet), 0, 0, 0, 0, 0, 0, -(E3*E3^(bet - 1)*MC^(((bet*thet - 1)*(thet - 1))/thet)*bet)/mcbar^(((bet*thet - 1)*(thet - 1))/thet), 0; 0, -(I*Q*gam)/(DA^(1/(alph - 1))*K), 0, 0, 0, 0, 0, 0, 0, 0, 0, Q*(gam*(ikbar - I/(DA^(1/(alph - 1))*K)) + 1), 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q, 0, 0, 0, 0, 0, 0, 0, 0, -E4; 0, 0, GDP, 0, (DA^(alph/(alph - 1))*H*K^alph*(alph - 1))/H^alph, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -C, -I, GDP, -GOV, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, R, -IR/E3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, IR/E3, 0; 0, 0, phiy*(rhoi - 1), 0, 0, 0, 0, 0, 1, phipi*(rhoi - 1), 0, 0, 0, phidy*(rhoi - 1), 0, 0, 0, 0, 0, 0, 0; 0, -I, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; C/(C - CL*DA^(1/(alph - 1))*kapp)^2 + (C*bet*kapp^2)/(C*kapp - C_p/DA_p^(1/(alph - 1)))^2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E1, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, E3, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (DA_p^(1/(alph - 1))*LAM_p*bet*(REN_p - Q_p*(delt + (gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p))^2)/2 - I_p/(DA_p^(1/(alph - 1))*K_p) + (I_p*(gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)) + 1))/(DA_p^(1/(alph - 1))*K_p) - 1)))/LAM, 0, 0, 0, 0, 0, 0, 0, 0, 0, E4; -C, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -IR, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -GDP, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -I, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -H, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -GDP/(DA^(1/(alph - 1))*YL), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, DY, 0, 0, 0, 0, 0, 0, 0; -C/(CL*DA^(1/(alph - 1))), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, DC, 0, 0, 0, 0, 0, 0; 0, -I/(DA^(1/(alph - 1))*IL), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, DI, 0, 0, 0, 0, 0; 0, 0, 0, 0, -H/HL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, DN, 0, 0, 0, 0; 0, 0, -GOV/GDP, GOV/GDP, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
fxp = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; K_p, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, (C_p*DA_p*bet*kapp)/(DA_p^(1/(alph - 1) + 1)*(C*kapp - C_p/DA_p^(1/(alph - 1)))^2*(alph - 1)); 0, 0, 0, 0, 0, 0, 0, 0, 0, -(DA_p*DA_p^(1/(alph - 1) - 1)*LAM_p*bet)/(alph - 1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; (DA_p^(1/(alph - 1))*LAM_p*Q_p*bet*(I_p/(DA_p^(1/(alph - 1))*K_p) + (I_p^2*gam)/(DA_p^(2/(alph - 1))*K_p^2) - (I_p*(gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)) + 1))/(DA_p^(1/(alph - 1))*K_p) + (I_p*gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)))/(DA_p^(1/(alph - 1))*K_p)))/LAM, 0, 0, 0, 0, 0, 0, 0, 0, (DA_p^(1/(alph - 1))*LAM_p*Q_p*bet*((DA_p*I_p)/(DA_p^(1/(alph - 1) + 1)*K_p*(alph - 1)) - (DA_p*I_p*(gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)) + 1))/(DA_p^(1/(alph - 1) + 1)*K_p*(alph - 1)) + (DA_p*I_p^2*gam)/(DA_p^(1/(alph - 1))*DA_p^(1/(alph - 1) + 1)*K_p^2*(alph - 1)) + (DA_p*I_p*gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)))/(DA_p^(1/(alph - 1) + 1)*K_p*(alph - 1))))/LAM - (DA_p*DA_p^(1/(alph - 1) - 1)*LAM_p*bet*(REN_p - Q_p*(delt + (gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p))^2)/2 - I_p/(DA_p^(1/(alph - 1))*K_p) + (I_p*(gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)) + 1))/(DA_p^(1/(alph - 1))*K_p) - 1)))/(LAM*(alph - 1)); 0, 0, CL_p, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, IRL_p, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, YL_p, 0, 0, 0, 0; 0, 0, 0, IL_p, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, HL_p, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
fyp = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -(C_p*bet*kapp)/(DA_p^(1/(alph - 1))*(C*kapp - C_p/DA_p^(1/(alph - 1)))^2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -DA_p^(1/(alph - 1))*LAM_p*bet, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -PI_p, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -(DA_p^(1/(alph - 1))*LAM_p*Q_p*bet*(I_p/(DA_p^(1/(alph - 1))*K_p) + (I_p^2*gam)/(DA_p^(2/(alph - 1))*K_p^2) - (I_p*(gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)) + 1))/(DA_p^(1/(alph - 1))*K_p) + (I_p*gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)))/(DA_p^(1/(alph - 1))*K_p)))/LAM, 0, 0, 0, 0, -(DA_p^(1/(alph - 1))*LAM_p*REN_p*bet)/LAM, 0, 0, 0, -(DA_p^(1/(alph - 1))*LAM_p*bet*(REN_p - Q_p*(delt + (gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p))^2)/2 - I_p/(DA_p^(1/(alph - 1))*K_p) + (I_p*(gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)) + 1))/(DA_p^(1/(alph - 1))*K_p) - 1)))/LAM, (DA_p^(1/(alph - 1))*LAM_p*Q_p*bet*(delt + (gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p))^2)/2 - I_p/(DA_p^(1/(alph - 1))*K_p) + (I_p*(gam*(ikbar - I_p/(DA_p^(1/(alph - 1))*K_p)) + 1))/(DA_p^(1/(alph - 1))*K_p) - 1))/LAM, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

eta = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 1; 0, 0; sigga, 0];
R = [];
