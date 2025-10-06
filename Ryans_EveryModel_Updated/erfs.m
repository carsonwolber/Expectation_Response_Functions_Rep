%**************************************************************
% SOLVE_LINEAR: Does some stuff.
%**************************************************************
addpath('../DSGE_tools');

nper  = 20;
nerf  = 16;

vnms = {'GDP', 'INV', 'CON', 'HRS'};


%%
%**************************************************************
% BLL ERF - TFP
%**************************************************************
clear *idx
[bll_param,bll_set] =bll_parameters;
load bll_obj

% %To get almost to RBC outcome...check why not more exact!
[f, fx, fy, fxp, fyp, eta]=bll_prog(struct2array(bll_param),struct2array(bll_set));

yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:da_idx-1];
xkidx = da_idx;
eqidx = 1:neq-1;
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
tfp_erf = yerf([dy_idx,di_idx,dc_idx,n_idx],:);
tfp_erf = tfp_erf./tfp_erf(1,1);

%%
%**************************************************************
% BLL ERF - GOV
%**************************************************************
yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:g_idx-1,g_idx+1:da_idx];
xkidx = g_idx;
eqidx = [1:neq-2,neq];
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
g_erf = yerf([dy_idx,di_idx,dc_idx,n_idx],:);
g_erf = g_erf./g_erf(1,1);

%%
%**************************************************************
% BLL ERF - MONEY
%**************************************************************
yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:q_idx-1,q_idx+1:da_idx];
xkidx = q_idx;
eqidx = [1:neq-3,neq-1:neq];
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
mon_erf = yerf([dy_idx,di_idx,dc_idx,n_idx],:);
mon_erf = mon_erf./mon_erf(1,1);

%%
%**************************************************************
% BLL ERF - EPSP
%**************************************************************
yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:mpt_idx-1,mpt_idx+1:da_idx];
xkidx = mpt_idx;
eqidx = [1:neq-4,neq-2:neq];
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
mpt_erf = yerf([dy_idx,di_idx,dc_idx,n_idx],:);
mpt_erf = mpt_erf./mpt_erf(1,1);

%%
%**************************************************************
% BLL ERF - EPSW
%**************************************************************
yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:mwt_idx-1,mwt_idx+1:da_idx];
xkidx = mwt_idx;
eqidx = [1:neq-5,neq-3:neq];
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
mwt_erf = yerf([dy_idx,di_idx,dc_idx,n_idx],:);
mwt_erf = mwt_erf./mwt_erf(1,1);
%%
%**************************************************************
% BLL ERF - ISP
%**************************************************************
yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:d_idx-1,d_idx+1:da_idx];
xkidx = d_idx;
eqidx = [1:neq-6,neq-4:neq];
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
isp_erf = yerf([dy_idx,di_idx,dc_idx,n_idx],:);
isp_erf = isp_erf./isp_erf(1,1);




%% Make Figure
f = figure;

p = cell(4,6);
for jj = 1:4
   s = subplot(2,2,jj);
   hold on 
   
   p{jj,1} = plot(0:nerf,g_erf(jj,:)  , '-r','linewidth',2);
   p{jj,2} = plot(0:nerf,isp_erf(jj,:), '-b','linewidth',2);
   p{jj,3} = plot(0:nerf,tfp_erf(jj,:), '-k','linewidth',2);
   p{jj,4} = plot(0:nerf,mon_erf(jj,:)  , '-g','linewidth',2);
   p{jj,5} = plot(0:nerf,mpt_erf(jj,:)  , '-','linewidth',2 , 'color', [.6 .4 .8]);
   p{jj,6} = plot(0:nerf,.2*mwt_erf(jj,:)  , '-','linewidth',2 , 'color', [.9 .8 .4]);
   
   
   s.XLim = [0,length(g_erf)-1];
   title(vnms{jj})
   if jj == 1
       xlabel('horizon')
   end
end

for jj = 1:4
    s = subplot(2,2,jj);
    plot(0:nerf, zeros(1,nerf+1), ':k');
    s.YLim = s.YLim;
end
legend('G','ISP','TFP', 'IR', 'PM', 'WM')

%% Now make it a layered figure
for ss = 2:6
    for jj = 1:4
        p{jj,ss}.Visible = false;
    end
end

saveas(f, '../../../slides/BC_macro_lunch/Figures/shock_compare1.eps', 'epsc');

for ss = 2:6
   
    
    
    for jj = 1:4
        p{jj,ss-1}.Color = .75+min(.2*p{jj,ss-1}.Color,1); 
        p{jj,ss}.Visible = true;
    end
    
    saveas(f, ['../../../slides/BC_macro_lunch/Figures/shock_compare' num2str(ss) '.eps'], 'epsc');

end