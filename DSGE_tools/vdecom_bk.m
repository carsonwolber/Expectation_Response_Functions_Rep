function [Vy, dOUT, mod_mom_blck,gx,hx,eta, out_pure,Ceq,dev]  =vdecom_bk(gx,hx,eta,bk_wghts,trend,var_idx,var_names,varargin)%(param,set,bk_per,bk_wghts,trend,var_idx,var_names,varargin)

bk_per = (length(bk_wghts)-1)/2;


if nargin>7
    prefs = varargin{1};
else
    prefs.ndec = 2;
end

sformat = ['%1.' num2str(prefs.ndec), 'f\t'];


ny = size(gx,1);
nx = size(hx,2);
neps = size(eta,2);



%Number of target variables in target moments.
S = selector(var_idx,ny+nx);
nvar = size(S,1);  
fmat = make_filter(nvar,bk_wghts,trend,0,1);


%Get sparse versions of eq. matrixes, setting near zeros to hard zeros.
hx    = nsparse(hx);
gx    = nsparse(gx);
gx_sm = nsparse(S*[gx;eye(nx)]);


%**************************************************************************
% THIS PART OF THE CODE COMPUTE THE MAIN MOMENT MATCHING LOSS FUNCTION
%**************************************************************************

mod_mom = zeros(neps,nvar);
for ll = 1:neps
    
    eta_tmp = zeros(size(eta));
    eta_tmp(:,ll) = eta(:,ll);
   % eta_tmp = eta;
    %Compute the covariance matrix sigX
    [~,sigX] = mom(gx_sm,hx,eta_tmp*eta_tmp',0);
    
    %Compue auto-covariance and orthogonality violation using recursion.
    sigXJ = sigX;
    mom_block = zeros(nvar,nvar*(2*bk_per+1));
    for kk = 0:(2*bk_per+1)
        
        %Cross moments
        mom_block(:,kk*nvar+(1:nvar)) = gx_sm*sigXJ*gx_sm';
        
        sigXJ = hx*sigXJ;
    end
    
    %Compute pure loss
    mod_mom_blck(:,:,ll) = reshape(fmat*mom_block(:),[nvar,nvar]);
    
    mod_mom(ll,:) = diag(mod_mom_blck(1:nvar,1:nvar,ll));

end


Vy = mod_mom./(repmat(sum(mod_mom),[size(eta,2),1]));


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

