% MODEL_SS - Return the steady state of the model computed analytically
%
% usage:
% 
% [ss, parameters] =model_ss(param)


function [ss,param,set] = rbc_ss(param,set)

%Upack parameters object
param_unpack


%BEGIN_EXTRACT_HERE




%Y and X vectors with SS values
Yss = zeros(1,100);
Xss = zeros(1,100);

ss = [Yss Xss];

%END_EXTRACT_HERE
