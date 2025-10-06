% VDECOM_TABLE - Make a variance decomp. table, print to screen.
%
% usage
%
% out = vedecom_table(gx,hx,eta,var_idx,var_names,prefs)


function Vy = vdecom_table(gx,hx,eta,var_idx,var_names,varargin)


if nargin>5
    prefs = varargin{1};
else
    prefs.ndec = 8;
    prefs.toscr = true;
end

sformat = ['%1.' num2str(prefs.ndec), 'f\t'];





[Vy,Vx] = variance_decomposition(gx,hx,eta);
neps = size(Vy,1);

Vy = [Vy,Vx];

Vy = Vy(:,var_idx);

if ~isfield(prefs, 'toscr') || prefs.toscr==true
    disp('Variance Decomposition')
    yidx = var_idx(1);
    tit_str = [];
    for j = 1:length(var_idx)
        idx = var_idx(j);
        tit_str = [tit_str, var_names{j}, '\t'];
    end
    
    
    disp(sprintf(tit_str))
    for j = 1:neps
        disp(sprintf(sformat, Vy(j,:)))
    end
end