% FEVD_TABLE - Make a variance decomp. table, print to screen.
%
% usage
%
% out = vedecom_table(gx,hx,eta,h,var_idx,var_name,prefs)


function Vy = fevd_table(gx,hx,eta,h,var_idx,var_name,varargin)


if nargin>6
    prefs = varargin{1};
else
    prefs.ndec = 2;
end

sformat = ['%1.' num2str(prefs.ndec), 'f\t'];




disp(['FEV Decomposition: ' var_name])
[Vy,Vx] = fevd(gx,hx,eta,h);
neps = size(Vy,1);


tit_str = [];
for j = 1:length(h)
    tit_str = [tit_str,'H=' num2str(h(j)) '\t'];
end


disp(sprintf(tit_str))
for j = 1:neps
    disp(sprintf(sformat, Vy(j,var_idx,:))) 
end