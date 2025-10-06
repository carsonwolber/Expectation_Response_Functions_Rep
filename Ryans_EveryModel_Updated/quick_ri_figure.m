%**************************************************************
% BLL ERF
%**************************************************************
clear *idx
[bll_param,bll_set] =bll_parameters;
load bll_obj

nerf = 500;

phipi = linspace(1.01,5,10);

sig_pi_ri = zeros(1,length(phipi));
sig_pi_ii = zeros(1,length(phipi));
for jj = 1:length(phipi)
    
    %**************************************************************
    % NK ERF
    %**************************************************************
    clear *idx
    [nk_param,nk_set] =nk_parameters;
    load nk_obj
    
    % Starndard parameters
    nk_set.phipi =  phipi(jj);       %taylor inflation
    nk_set.phiy  =  0.5;       %taylor output level
    nk_set.phidy =  0  ;       %taylor rule on output growth
    nk_set.rhoi  =  0.5;       %taylor smoothing
    
    
    [f, fx, fy, fxp, fyp, eta]=nk_prog(struct2array(nk_param),struct2array(nk_set));

    [gx,hx] = gx_hx_alt(fy,fx,fyp,fxp);
    
    sigY = mom(gx,hx,eta*eta');
    
    sig_pi_ii(jj) = sigY(pi_idx,pi_idx);
    
    
    yidx = [gdp_idx:e2_idx];
    xidx = [irl_idx:da_idx-1];
    xkidx = da_idx;
    eqidx = 1:neq-1;
    yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
    

    sig_pi_ri(jj) = sum(yerf(pi_idx,:).^2);
    
    
end