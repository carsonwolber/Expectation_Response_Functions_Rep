% PARAMETETRS - This function returns a parameter structure to use in the model solution.


function [param,set] = nk_parameters()
                                                       
set.adiff      = 0;
set.approx_deg = 1;


%parameters B&S calibrate
set.bet = .99;            %discount factor
set.alph = 0.1859;        %capital share
set.thet = .8770;%.8770;  %Price change prob
mu_p = .3;
set.etap = (1+mu_p)/mu_p;   %Output elasticity
set.zeta = 2.0871;   %Labor supply elasiticity

%parameters B&S estimate
set.rhoe  = 0.94;         %??
param.rhoga = 0.73;         %autocorrelation of growth rate
param.sigga = 0.17;         %shock to growth rate
param.sigea = 0.58;         %shock to level
param.siges = 0.13;         %noise in signal

% Starndard parameters
 set.phipi =  1.5;       %taylor inflation
 set.phiy  =  0.5;       %taylor output level
 set.phidy =  0  ;       %taylor rule on output growth
 set.rhoi  =  0.5;       %taylor smoothing
set.gstar = 1;
