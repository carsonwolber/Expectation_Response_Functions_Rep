% PARAMETETRS - This function returns a parameter structure to use in the model solution.


function [param,set] = rbc_parameters()
                                                       
set.adiff      = 0;
set.approx_deg = 1;


%Parameters of physical environment
set.bet  = .985;          %discount factor
set.alph = 0.36;          %capital share
set.del  = 0.0125;
set.chi  = 1;
set.gam  = exp(0);

%Parameters of exog process
param.rho  = 0;
param.siga = 1;