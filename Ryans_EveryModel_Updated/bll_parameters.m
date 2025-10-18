% PARAMETETRS - This function returns a parameter structure to use in the model solution.
function [param,set] = bll_parameters()                                                
set.adiff      = 0;
set.approx_deg = 1;


%***************************************
% PARAMETERS FOR ECONOMY
%***************************************
%Parameters BLL info
param.rho = .9426;
param.sigtmp = 1.1977;
param.sigv = 1.4738;

% param.rhox = .9426;
% param.rhoz = .9426;
% param.sige = param.sigtmp*(1-param.rho);
% param.sign = sqrt(param.rho)*param.sigtmp;
% param.sigv = 1.4738;


%economic parameters BLL, calibrated?
set.bet = .99;            %discount factor
set.delt = .025;         %Depreciation
set.gam   = 1.00;         %Growth Rate
set.psii  = .21875;           %Goverment share
set.muw = .05;  

%Estimatied
param.h     = 0.5262;          %habit
param.alph  = 0.1859;        %labor share
param.zeta  = 2.0871;          %frisch
param.xi   = 3.4919;         %capacity cost
param.chi = 4.3311;          %adjustment cost
param.thet = 0.8770;          %calvo price
param.thetw = 0.8690;         %calvo wage
      

%Other estimated exogenous processes
param.rhod = .4641;
param.sigd = 11.098;
param.rhop = .7722;
param.phip = -.4953;
param.sigp = .1778;
param.rhow = .9530;
param.phiw = -.9683;
param.sigw = .3057;


param.rhor = 0.5583;
param.rhoq = 0.0413;
param.sigq = 0.3500;
param.rhog = 0.9972;
param.sigg = 0.2877;

% Starndard parameters
%  set.phipi =  1.5;       %taylor inflation
%  set.phiy  =  0.5;       %taylor output level
%  set.phidy =  0  ;       %taylor rule on output growth
%  set.rhoi  =  0.5;       %taylor smoothing

%BS parameters
% set.phipi = 1.31;         %taylor rule on inflation
% set.phiy  = 0   ;         %taylor rule on output level
% set.phidy = 0.94;         %taylor rule on output growth
% set.rhoi  = 0.66;         %taylor smoothing

%BLL parameters
 set.phipi =  1.0137;       %taylor inflation
 set.phiy  = 0.0050;        %taylor output level
 set.phidy = 0;             %taylor rule on output growth
 set.rhoi  = 0.5583;        %taylor smoothing




%Determined in equilbrium
set.cybar  = NaN;      
set.iybar = NaN;
set.rkpybar = NaN;
set.kap = NaN;
set.kapw = NaN;
