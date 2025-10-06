function [f fx fy fxp fyp eta R set dgam_dtheta deta_dtheta dR_dtheta xlag ylag] = model_prog(param, set)

%Assign parameter values to named variables.
rho = param(1);
siga = param(2);

%Assign set values to named variables.
adiff = set(1);
approx_deg = set(2);
bet = set(3);
alph = set(4);
del = set(5);
chi = set(6);
gam = set(7);

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
%Compute Steady State
GAM= Xss(1);
A= Xss(2);
K= Xss(3);
CL= Xss(4);
IL= Xss(5);
YL= Xss(6);
GDP= Yss(1);
C= Yss(2);
H= Yss(3);
W= Yss(4);
R= Yss(5);
I= Yss(6);
DC= Yss(7);
DI= Yss(8);
DY= Yss(9);
GAM_p= Xss(1);
A_p= Xss(2);
K_p= Xss(3);
CL_p= Xss(4);
IL_p= Xss(5);
YL_p= Xss(6);
GDP_p= Yss(1);
C_p= Yss(2);
H_p= Yss(3);
W_p= Yss(4);
R_p= Yss(5);
I_p= Yss(6);
DC_p= Yss(7);
DI_p= Yss(8);
DY_p= Yss(9);

%Evaluate F.
f = [C - W/chi, 1 - (C*GAM_p^(1/(alph - 1))*bet*(R_p - del + 1))/C_p, R - A*GAM*alph*(K/H)^(alph - 1), W + A*GAM^(alph/(alph - 1))*(K/H)^alph*(alph - 1), I/GAM^(1/(alph - 1)) - K*(del - 1) - K_p/GAM^(1/(alph - 1)), C + I - A*GAM^(alph/(alph - 1))*H^(1 - alph)*K^alph, GDP - A*GAM^(alph/(alph - 1))*H^(1 - alph)*K^alph, CL_p - C, IL_p - I, YL_p - GDP, DC - C/(CL*GAM^(1/(alph - 1))), DI - I/(GAM^(1/(alph - 1))*IL), DY - GDP/(GAM^(1/(alph - 1))*YL), log(GAM_p/gam) - rho*log(GAM/gam), log(A_p) - (9*log(A))/10];
%Evaluate derivative expressions.
fx = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -A*GAM*alph*(K/H)^(alph - 1), -A*GAM*alph*(K/H)^(alph - 1), -(A*GAM*K*alph*(K/H)^(alph - 2)*(alph - 1))/H, 0, 0, 0; A*GAM*GAM^(alph/(alph - 1) - 1)*alph*(K/H)^alph, A*GAM^(alph/(alph - 1))*(K/H)^alph*(alph - 1), (A*GAM^(alph/(alph - 1))*K*alph*(K/H)^(alph - 1)*(alph - 1))/H, 0, 0, 0; (GAM*K_p)/(GAM^(1/(alph - 1) + 1)*(alph - 1)) - (GAM*I)/(GAM^(1/(alph - 1) + 1)*(alph - 1)), 0, -K*(del - 1), 0, 0, 0; -(A*GAM*GAM^(alph/(alph - 1) - 1)*H^(1 - alph)*K^alph*alph)/(alph - 1), -A*GAM^(alph/(alph - 1))*H^(1 - alph)*K^alph, -A*GAM^(alph/(alph - 1))*H^(1 - alph)*K*K^(alph - 1)*alph, 0, 0, 0; -(A*GAM*GAM^(alph/(alph - 1) - 1)*H^(1 - alph)*K^alph*alph)/(alph - 1), -A*GAM^(alph/(alph - 1))*H^(1 - alph)*K^alph, -A*GAM^(alph/(alph - 1))*H^(1 - alph)*K*K^(alph - 1)*alph, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; (C*GAM)/(CL*GAM^(1/(alph - 1) + 1)*(alph - 1)), 0, 0, C/(CL*GAM^(1/(alph - 1))), 0, 0; (GAM*I)/(GAM^(1/(alph - 1) + 1)*IL*(alph - 1)), 0, 0, 0, I/(GAM^(1/(alph - 1))*IL), 0; (GAM*GDP)/(GAM^(1/(alph - 1) + 1)*YL*(alph - 1)), 0, 0, 0, 0, GDP/(GAM^(1/(alph - 1))*YL); -rho, 0, 0, 0, 0, 0; 0, -9/10, 0, 0, 0, 0];
fy = [0, C, 0, -W/chi, 0, 0, 0, 0, 0; 0, -(C*GAM_p^(1/(alph - 1))*bet*(R_p - del + 1))/C_p, 0, 0, 0, 0, 0, 0, 0; 0, 0, (A*GAM*K*alph*(K/H)^(alph - 2)*(alph - 1))/H, 0, R, 0, 0, 0, 0; 0, 0, -(A*GAM^(alph/(alph - 1))*K*alph*(K/H)^(alph - 1)*(alph - 1))/H, W, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, I/GAM^(1/(alph - 1)), 0, 0, 0; 0, C, (A*GAM^(alph/(alph - 1))*H*K^alph*(alph - 1))/H^alph, 0, 0, I, 0, 0, 0; GDP, 0, (A*GAM^(alph/(alph - 1))*H*K^alph*(alph - 1))/H^alph, 0, 0, 0, 0, 0, 0; 0, -C, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -I, 0, 0, 0; -GDP, 0, 0, 0, 0, 0, 0, 0, 0; 0, -C/(CL*GAM^(1/(alph - 1))), 0, 0, 0, 0, DC, 0, 0; 0, 0, 0, 0, 0, -I/(GAM^(1/(alph - 1))*IL), 0, DI, 0; -GDP/(GAM^(1/(alph - 1))*YL), 0, 0, 0, 0, 0, 0, 0, DY; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0];
fxp = [0, 0, 0, 0, 0, 0; -(C*GAM_p*GAM_p^(1/(alph - 1) - 1)*bet*(R_p - del + 1))/(C_p*(alph - 1)), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, -K_p/GAM^(1/(alph - 1)), 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, CL_p, 0, 0; 0, 0, 0, 0, IL_p, 0; 0, 0, 0, 0, 0, YL_p; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0];
fyp = [0, 0, 0, 0, 0, 0, 0, 0, 0; 0, (C*GAM_p^(1/(alph - 1))*bet*(R_p - del + 1))/C_p, 0, 0, -(C*GAM_p^(1/(alph - 1))*R_p*bet)/C_p, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0];

eta = [siga; 0; 0; 0; 0; 0];
R = [];
