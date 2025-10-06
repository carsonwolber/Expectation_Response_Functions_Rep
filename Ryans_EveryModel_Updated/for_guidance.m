%**************************************************************
% SOLVE_LINEAR: Does some stuff.
%**************************************************************
addpath('../DSGE_tools');

nerf  = 24;

vnms = {'GDP', 'CON', 'PI', 'IR'};

%%
%**************************************************************
% NK ERF
%**************************************************************
clear *idx
[nk_param,nk_set] =nk_parameters;
load nk_obj
[f, fx, fy, fxp, fyp, eta]=nk_prog(struct2array(nk_param),struct2array(nk_set));

yidx = [gdp_idx:e2_idx];
xidx = [irl_idx:ei_idx-1,ei_idx+1:da_idx];
xkidx = ei_idx;
eqidx = [1:neq-3,neq-1:neq];
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
nk_erf = yerf([dy_idx,dy_idx,pi_idx,ir_idx],:);

%%
%**************************************************************
% BS ERF
%**************************************************************
clear *idx
[bs_param,bs_set] =bs_parameters;
load bs_obj

% %To get almost to RBC outcome...check why not more exact!
% bs_set.gam = .0001;
% bs_set.etta = 100000;
% bs_set.thet = .001;
% bs_set.kap = .0001;
% bs_set.xi = 10000;
% bs_set.Gshr = .00001;
[f, fx, fy, fxp, fyp, eta]=bs_prog(struct2array(bs_param),struct2array(bs_set));

yidx = [c_idx:e4_idx];
xidx = [k_idx:ei_idx-1,ei_idx+1:da_idx];
xkidx = ei_idx;
eqidx = [1:neq-3,neq-1:neq];
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
bs_erf = yerf([dy_idx,dc_idx,pi_idx,ir_idx],:);

%%
%**************************************************************
% BLL ERF
%**************************************************************
clear *idx
[bll_param,bll_set] =bll_parameters;
load bll_obj

% %To get almost to RBC outcome...check why not more exact!
[f, fx, fy, fxp, fyp, eta]=bll_prog(struct2array(bll_param),struct2array(bll_set));
yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:q_idx-1,q_idx+1:da_idx];
xkidx = q_idx;
eqidx = [1:neq-3,neq-1:neq];
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
bll_erf = yerf([dy_idx,dc_idx,pi_idx,irt_idx],:);


%% Make Figure
f = figure;
for jj = 1:4
   s = subplot(2,2,jj);
   plot(0:nerf,nk_erf(jj,:), '-k','linewidth',2);
   hold on 
   plot(0:nerf,bs_erf(jj,:), '-g','linewidth',2);
   plot(0:nerf,bll_erf(jj,:), '-r','linewidth',2);
   
   %plot(0:length(rbc_erf)-1,bs_erf(jj,:), '-r','linewidth',2);
   
   
   s.XLim = [0,nerf];
   title(vnms{jj})
   if jj == 1
       xlabel('horizon')
   end
end
legend('NK','BS','MS')

saveas(f, '../../../slides/BC_macro_lunch/Figures/fg_compare.eps', 'epsc');