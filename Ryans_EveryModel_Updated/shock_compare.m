%**************************************************************
% SOLVE_LINEAR: Does some stuff.
%**************************************************************
addpath('DSGE_tools');
nper  = 20;
nerf  = 5000;

vnms = {'GDP',  'CON','INV', 'HRS'};


%%
%**************************************************************
% BLL ERF - TFP
%**************************************************************
clear *idx
[bll_param,bll_set] =bll_parameters;
load bll_obj

erf_idx = [dy_idx,dc_idx,pi_idx,irt_idx];

% Starndard parameters
bll_set.phipi =  1.5;       %taylor inflation
bll_set.phiy  =  0.5;       %taylor output level
bll_set.phidy =  0  ;       %taylor rule on output growth
bll_set.rhoi  =  0.5;       %taylor smoothing

yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:al_idx];
eqidx = 1:length([yidx,xidx]);

% %To get almost to RBC outcome...check why not more exact!
[f, fx, fy, fxp, fyp, eta]=bll_prog(struct2array(bll_param),struct2array(bll_set));
%%
%**************************************************************
% BLL ERF - GOV
%**************************************************************
xkidx = g_idx;
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
g_erf = yerf(erf_idx,:);
g_erf = g_erf./g_erf(1,1);

[all_sfw,all_mhe,all_lpf32,all_lpf100] = erf_stats(g_erf(1,:),GAM,ZZ(dy_idx,:),TQA,TQB);


%%
%**************************************************************
% BLL ERF - TFP
%**************************************************************
xkidx = a_idx;
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
tfp_erf = yerf(erf_idx,:);
tfp_erf = tfp_erf./tfp_erf(1,1);

[all_sfw(2),all_mhe(2),all_lpf32(2),all_lpf100(2)] = erf_stats(tfp_erf(1,:),GAM,ZZ(dy_idx,:),TQA,TQB);


%**************************************************************
% BLL ERF - ISP
%**************************************************************
xkidx = d_idx;
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
isp_erf = yerf(erf_idx,:);
isp_erf = isp_erf./isp_erf(1,1);
[all_sfw(3),all_mhe(3),all_lpf32(3),all_lpf100(3)] = erf_stats(isp_erf(1,:),GAM,ZZ(dy_idx,:),TQA,TQB);



%%
%**************************************************************
% BLL ERF - MONEY
%**************************************************************
xkidx = q_idx;
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
mon_erf = yerf(erf_idx,:);
mon_erf = mon_erf./mon_erf(1,1);

[all_sfw(4),all_mhe(4),all_lpf32(4),all_lpf100(4)] = erf_stats(mon_erf(1,:),GAM,ZZ(dy_idx,:),TQA,TQB);


%%
%**************************************************************
% BLL ERF - EPSP
%**************************************************************
xkidx = mpt_idx;
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
mpt_erf = yerf(erf_idx,:);
mpt_erf = mpt_erf./mpt_erf(1,1);
[all_sfw(5),all_mhe(5),all_lpf32(5),all_lpf100(5)] = erf_stats(mpt_erf(1,:),GAM,ZZ(dy_idx,:),TQA,TQB);


%%
%**************************************************************
% BLL ERF - EPSW
%**************************************************************
xkidx = mwt_idx;
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
mwt_erf = yerf(erf_idx,:);
mwt_erf = mwt_erf./mwt_erf(1,1);

[all_sfw(6),all_mhe(6),all_lpf32(6),all_lpf100(6)] = erf_stats(mwt_erf(1,:),GAM,ZZ(dy_idx,:),TQA,TQB);

%% **************************************************************
% BLL ERF - EPSW
%**************************************************************
shock_list = [g_idx,gam_idx,d_idx,q_idx,mpt_idx,mwt_idx];

nphi = 25;
phi_grid = linspace(1.05,4,nphi);
ystat = zeros(4,nphi,6); ystat2 = 0;
bll_set.phiy  = 0.5;        %taylor output level
bll_set.phidy = 0;             %taylor rule on output growth
bll_set.rhoi  = 0.5;        %taylor smoothing

for jj =1:length(phi_grid)
    bll_set.phipi = phi_grid(jj);
    [f, fx, fy, fxp, fyp, eta]=bll_prog(struct2array(bll_param),struct2array(bll_set));

    for ii = 1:6
        [yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,shock_list(ii),eqidx,nerf);

        [~,~,~,ystat(:,jj,ii)] = erf_stats(yerf(erf_idx,:),GAM,ZZ(erf_idx,:),TQA,TQB);
    end
end
   
%% Make Figure
f = figure;

p = cell(1,6);
ls = {'-r','-k','-b','-g','-','-'};
lcolor = {'r','k','b','g',[.6 .4 .8],[0    0.4471    0.7412]};

for jj = 1:4
    s = subplot(2,2,jj);
    for ii = 1:6

        hold on

        p{ii} = plot(phi_grid,ystat(jj,:,ii) , ls{ii}, 'linewidth', 2 , 'color', lcolor{ii});



        if jj == 1
            xlabel('\phi_\pi')
            ylabel('lpf(100)')
        end
    end
      title(vnms{jj})
end

legend('G','TFP','ISP', 'IR', 'PM', 'WM')

saveas(f, ['figures_tables/shock_compare_tight' num2str(ss) '.eps'], 'epsc');


exportgraphics(f, 'figures_tables/shock_compare_tight.jpg','Resolution',300);



%% Now make table commented due to lack of table format
% output_dat = zeros(6,4);
% 
% output_dat(:,1) = all_sfw(:);
% output_dat(:,2) = all_mhe(:);
% output_dat(:,3) = all_lpf32(:);
% output_dat(:,4) = all_lpf100(:);
% 
% 
% table_insert('../../../paper/figures_tables/main_table_all.text', '../../../paper/figures_tables/main_table_all.tex',...
%     output_dat, {'%0.2f','%0.2f','%0.2f'});