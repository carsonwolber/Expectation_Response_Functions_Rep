% PARAMETETRS - This function returns a parameter structure to use in the model solution.
function [param,set] = bs_parameters()
                                               
set.adiff      = 0;
set.approx_deg = 1;



%parameters B&S calibrate
set.bet = .99;            %discount factor
set.alph = 0.36;          %capital share
set.delt = .03;            %depreciation rate
set.gstar = exp(.0033);   %steady-state TFP growth
set.Gshr = .2;            %govm't share
set.rhog = 0.95;          %peristence of govenment rate
set.sigg = .25;           %govt shock.

%parameters B&S estimate
set.rhoe  = 0.94;         %??
param.rhoga = 0.73;         %autocorrelation of growth rate
param.sigga = 0.17;         %shock to growth rate
param.sigea = 0.58;         %shock to level
param.siges = 0.13;         %noise in signal

% param.rhoga = 0.735775021120853;         %autocorrelation of growth rate
% param.sigga = 0.621548429963160;         %shock to growth rate
% param.sigea = 0.769525358338735;         %shock to level
% param.siges = 0.011925268462234;         %noise in signal

% Starndard parameters
%  set.phipi =  1.5;       %taylor inflation
%  set.phiy  =  0.5;       %taylor output level
%  set.phidy =  0  ;       %taylor rule on output growth
%  set.rhoi  =  0.5;       %taylor smoothing

%BS parameters
set.phipi = 1.31;         %taylor rule on inflation
set.phiy  = 0;
set.phidy  = 0.94;         %taylor rule on output growth
set.rhoi  = 0.66;         %taylor smoothing

%BLL parameters
%  set.phipi =  1.0137;       %taylor inflation
%  set.phiy  = 0.0050;        %taylor output level
%  set.phidy = 0;             %taylor rule on output growth
%  set.rhoi  = 0.5583;        %taylor smoothing






set.kapp  = 0.31;          %habit
set.gam   = 0.16/(set.delt+set.gstar^(1/(1-set.alph))-1); %adjustment cost elasticity (denom aligns with B&S code)
set.etta  = 1.32;          %labor supply elasticity (don't confuse with eta, covariance matrix of structural shocks)
set.sigi  = .5*0.21;          %shocks to taylor rule?
set.xi    = 13.71;         %elatisticy of substitution in aggregator
set.thet  = 0.76;          %calvo prop of not changing


%Determined in equilbrium
set.ybar  = NaN;      
set.mcbar = NaN;
set.ikbar = NaN;
set.irbar = NaN;
set.gdp   = NaN;






