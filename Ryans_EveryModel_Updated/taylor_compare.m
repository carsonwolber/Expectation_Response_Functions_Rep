%**************************************************************
% SOLVE_LINEAR: Does some stuff.
%**************************************************************
addpath('DSGE_tools');

nper  = 20;
nerf  = 5000;

vnms = {'GDP', 'CON', 'INV', 'HRS'};



%%
%**************************************************************
% NK MODEL
%**************************************************************
clear *idx
[nk_param,nk_set] =nk_parameters;
load nk_obj

yidx = [gap_idx:dy_idx];
xidx = [irl_idx:a_idx-2];
xkidx = gam_idx;
eqidx = 1:length([yidx,xidx]);


% Starndard parameters
nk_set.phipi =  1.5;       %taylor inflation
nk_set.phiy  =  0.5;       %taylor output level
nk_set.phidy =  0  ;       %taylor rule on output growth
nk_set.rhoi  =  0.5;       %taylor smoothing

[f, fx, fy, fxp, fyp, eta]=nk_prog(struct2array(nk_param),struct2array(nk_set));
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
std_erf = yerf([dy_idx,dy_idx,dy_idx,gap_idx],:);
std_erf(3,:) = NaN; %No Investment
[std1_sfw,std1_mhe,std1_lpf32,std1_lpf100] = erf_stats(std_erf,GAM,ZZ([dy_idx,dy_idx,dy_idx,gap_idx],:),TQA,TQB);


%BLL parameters
nk_set.phipi =  1.0137;       %taylor inflation
nk_set.phiy  = 0.0050;        %taylor output level
nk_set.phidy = 0;             %taylor rule on output growth
nk_set.rhoi  = 0.5583;        %taylor smoothing


[f, fx, fy, fxp, fyp, eta]=nk_prog(struct2array(nk_param),struct2array(nk_set));
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
bll_erf = yerf([dy_idx,dy_idx,dy_idx,gap_idx],:);
bll_erf(3,:) = NaN; %No Investment
[bll1_sfw,bll1_mhe,bll1_lpf32,bll1_lpf100] = erf_stats(bll_erf,GAM,ZZ([dy_idx,dy_idx,dy_idx,gap_idx],:),TQA,TQB);



%%
%**************************************************************
% BLL MODEL
%**************************************************************
clear *idx
[bll_param,bll_set] =bll_parameters;
load bll_obj

yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:al_idx];
xkidx = gam_idx;
eqidx = 1:length([yidx,xidx]);

% Starndard parameters
bll_set.phipi =  1.5;       %taylor inflation
bll_set.phiy  =  0.5;       %taylor output level
bll_set.phidy =  0  ;       %taylor rule on output growth
bll_set.rhoi  =  0.5;       %taylor smoothing

% %To get almost to RBC outcome...check why not more exact!
[f, fx, fy, fxp, fyp, eta]=bll_prog(struct2array(bll_param),struct2array(bll_set));
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
std_erf(:,:,3) = yerf([dy_idx,dc_idx,di_idx,n_idx],:);
[std3_sfw,std3_mhe,std3_lpf32,std3_lpf100] = erf_stats(std_erf(:,:,3),GAM,ZZ([dy_idx,dc_idx,di_idx,n_idx],:),TQA,TQB);

%BLL parameters
bll_set.phipi =  1.0137;       %taylor inflation
bll_set.phiy  = 0.0050;        %taylor output level
bll_set.phidy = 0;             %taylor rule on output growth
bll_set.rhoi  = 0.5583;        %taylor smoothing

[f, fx, fy, fxp, fyp, eta]=bll_prog(struct2array(bll_param),struct2array(bll_set));
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
bll_erf(:,:,3) = yerf([dy_idx,dc_idx,di_idx,n_idx],:);
[bll3_sfw,bll3_mhe,bll3_lpf32,bll3_lpf100] = erf_stats(bll_erf(:,:,3),GAM,ZZ([dy_idx,dc_idx,di_idx,n_idx],:),TQA,TQB);





%% Make Figure
f = figure;

p = cell(4,6);
for jj = 1:4
   s = subplot(2,2,jj);
    hold on 
   p{jj,1} = plot(0:jplot-1,std_erf(jj,1:jplot,1), '-b','linewidth',2, 'color',[    0.0275    0.6706    0.5843]);
   p{jj,3} = plot(0:jplot-1,std_erf(jj,1:jplot,3), '-','linewidth',2, 'color',[  0.5333    0.0314    0.6118]);
   p{jj,4} = plot(0:jplot-1,bll_erf(jj,1:jplot,1), '-xb','linewidth',2, 'color',[    0.0275    0.6706    0.5843]);
   p{jj,6} = plot(0:jplot-1,bll_erf(jj,1:jplot,3), '-x','linewidth',2, 'color',[  0.5333    0.0314    0.6118]);
   
   
   
   s.XLim = [0,jplot-1];
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
%legend('NK', 'NK-RBC', 'SW', 'NK''', 'NK-RBC''', 'SW''')

legend('New Keyn.: stand. policy',  'Med. Scale: stand. policy', 'New Keyn.: accom. policy',  'Med. Scale: accom. policy', 'Location','SouthEast')

%% Now make it a layered figure
for ss = 4:6
    for jj = 1:4
        p{jj,ss}.Visible = 'off';
    end
end

saveas(f, 'figures_tables/taylor_compare1.eps', 'epsc');

for ss = [1,3]
    for jj = 1:4
        p{jj,ss}.Color = .75+min(.2*p{jj,ss}.Color,1); 
        p{jj,ss+3}.Visible = 'on';
    end
end
saveas(f, 'figures_tables/taylor_compare2.eps', 'epsc');


exportgraphics(f, 'figures_tables/taylor_compare2.jpg','Resolution',300);


%% Make table – Commented out due to current lack of table template


%
% tabidx = 1;  %Output
% %tabidx = 2;  %Consumption   
% %tabidx = 3;  %Investment
% %tabidx = 4;  %hours
% 
% output_dat = zeros(6,4);
% 
% output_dat(:,1) = [rbc_sfw(tabidx),jr_sfw(tabidx),std1_sfw(tabidx),std3_sfw(tabidx),bll1_sfw(tabidx),bll3_sfw(tabidx)]';
% output_dat(:,2) = [rbc_mhe(tabidx),jr_mhe(tabidx),std1_mhe(tabidx),std3_mhe(tabidx),bll1_mhe(tabidx),bll3_mhe(tabidx)]';
% output_dat(:,3) = [rbc_lpf32(tabidx),jr_lpf32(tabidx),std1_lpf32(tabidx),std3_lpf32(tabidx),bll1_lpf32(tabidx),bll3_lpf32(tabidx)]';
% output_dat(:,4) = [rbc_lpf100(tabidx),jr_lpf100(tabidx),std1_lpf100(tabidx),std3_lpf100(tabidx),bll1_lpf100(tabidx),bll3_lpf100(tabidx)]';
% 
% cons_dat = zeros(6,4);
% 
% cons_dat(:,1) = [rbc_sfw(2),jr_sfw(2),std1_sfw(2),std3_sfw(2),bll1_sfw(2),bll3_sfw(2)]';
% cons_dat(:,2) = [rbc_mhe(2),jr_mhe(2),std1_mhe(2),std3_mhe(2),bll1_mhe(2),bll3_mhe(2)]';
% cons_dat(:,3) = [rbc_lpf32(2),jr_lpf32(2),std1_lpf32(2),std3_lpf32(2),bll1_lpf32(2),bll3_lpf32(2)]';
% cons_dat(:,4) = [rbc_lpf100(2),jr_lpf100(2),std1_lpf100(2),std3_lpf100(2),bll1_lpf100(2),bll3_lpf100(2)]';
% 
% 
% table_insert('figures_tables/main_table.text', 'figures_tables/main_table.text',...
%     [output_dat,cons_dat], {'%0.2f','%0.2f','%0.2f'});
% 
% 
% 
% 
% table_insert('figures_tables/main_table_cons.text', 'figures_tables/main_table_cons.tex',...
%     cons_dat, {'%0.2f','%0.2f','%0.2f'});
