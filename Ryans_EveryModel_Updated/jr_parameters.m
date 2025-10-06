% PARAMETETRS - This function returns a parameter structure to use in the model solution.


function [param,set] = parameters()
                                                       
set.adiff      = 0;
set.approx_deg = 1;


%***************************************
% PARAMETERS FOR ECONOMY
%***************************************
param.sig     = 1.00000;   %Risk aversion
param.thet    = 1.4;     %Labor elasticity parameter
param.bet     = .985;  %.991 %Beta and R are separate things now
param.rworld  = 1/param.bet;
param.hbar     = 0.2;     %Steady-state labor
param.gamm     = .001;

param.alph    = 1-.64;    %Capital Share
param.phi     = 1.3;      %Investment adjustment cost parameter

param.delt0   = 0.0125;    %Depreciation rate  
param.delcurv = .15;
param.delt1   = NaN;
param.delt2   = NaN;

param.open = 0;
param.dshr = 0*.000001;
param.psii_ln = .01;

%Exogenous Parameters
param.rhox  = .0;     %AR coefficient of gamma_x
param.rhoa  = .60;    %AR coefficient of A
param.rhor  = .50;    %AR coefficient of RW

param.sigxn  = .01;     %News shock, permanent
param.sigxs  = .01;     %Suprise shock, permanent
param.siga   = .01;     %Tmp A shock
param.sigr   = .015;    %R shock



%Determined in equilbrium
param.psii  = NaN;      %Parameter of endogenous discount factor  
param.kbar  = NaN;      %Level of ss capital
param.dbar  = NaN;      %Level of ss debt
param.abar  = NaN;
param.rbar  = NaN;      %Actual domestic interest rate
param.dtarg = NaN;      %Share of debt holding that equates world and domstic interest rate
param.e2bar = NaN;
param.e3bar = NaN;
param.al = NaN;
param.gl = NaN;
param.rl = NaN;
param.ah = NaN;
param.gh = NaN;
param.rh = NaN;
