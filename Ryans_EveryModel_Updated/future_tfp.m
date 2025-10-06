%**************************************************************
% SOLVE_LINEAR: Does some stuff.
%**************************************************************
addpath('DSGE_tools/');

nper  = 20;
nerf  = 5000;

vnms = {'GDP', 'CON', 'INV', 'HRS'};

% Standard Solves
% [f, fx, fy, fxp, fyp, eta]=rbc_prog(struct2array(rbc_param),struct2array(rbc_set));
% [gx,hx]= gx_hx_alt(fy,fx,fyp,fxp);
% 
% ir1 = ir(gx,hx,eta,16)
% csum_idx = [dc_idx,di_idx,gam_idx];
% ir1(:,csum_idx) = cumsum(ir1(:,csum_idx));


%%
%**************************************************************
% RBC ERF
%**************************************************************
clear *idx
[rbc_param,rbc_set] = rbc_parameters;
load rbc_obj
%rbc_param.gam = 1;
yidx  = gdp_idx:dy_idx;
xidx  = gam_idx+2:yl_idx;
xkidx = gam_idx;
eqidx = 1:neq-2;


[f, fx, fy, fxp, fyp, eta]=rbc_prog(struct2array(rbc_param),struct2array(rbc_set));
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
rbc_erf = yerf([dy_idx,dc_idx,di_idx,h_idx],:);

[rbc_sfw,rbc_mhe,rbc_lpf32,rbc_lpf100] = erf_stats(rbc_erf,GAM,ZZ([dy_idx,dc_idx,di_idx,h_idx],:),TQA,TQB);


%%
%**************************************************************
% JR ERF
%**************************************************************
clear *idx
[jr_param,jr_set] =jr_parameters;
load jr_obj

yidx = [c_idx:r_idx];
xidx = [k_idx:gam_idx-1,gam_idx+2:d_idx];
xkidx = gam_idx;
eqidx = 1:neq-2;

[f, fx, fy, fxp, fyp, eta]=jr_prog(struct2array(jr_param),struct2array(jr_set));
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
jr_erf = yerf([dy_idx,dc_idx,di_idx,h_idx],:);
[jr_sfw,jr_mhe,jr_lpf32,jr_lpf100] = erf_stats(jr_erf,GAM,ZZ([dy_idx,dc_idx,di_idx,h_idx],:),TQA,TQB);

%%
%**************************************************************
% JR ERF
%**************************************************************
clear *idx
[jr_param,jr_set] =jr_parameters;
load jr_obj

jr_param.phi = 0.01;

yidx = [c_idx:r_idx];
xidx = [k_idx:gam_idx-1,gam_idx+2:d_idx];
xkidx = gam_idx;
eqidx = 1:neq-2;

[f, fx, fy, fxp, fyp, eta]=jr_prog(struct2array(jr_param),struct2array(jr_set));
[yerf,~,GAM,ZZ,TQA,TQB] = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
jr2_erf = yerf([dy_idx,dc_idx,di_idx,h_idx],:);
[jr2_sfw,jr2_mhe,jr2_lpf32,jr2_lpf100] = erf_stats(jr_erf,GAM,ZZ([dy_idx,dc_idx,di_idx,h_idx],:),TQA,TQB);






%% Make Figure
f = figure; f.PaperPosition = .5*[0.7500 3 7.0000 5];
f.Units = 'inches';

p = cell(4,5);
for jj = 1:4
   s = subplot(2,2,jj);
   p{jj,1} = plot(0:jplot-1,rbc_erf(jj,1:jplot), '-k','linewidth',2);
   hold on 
   p{jj,2} = plot(0:jplot-1,jr_erf(jj,1:jplot), '-','linewidth',2, 'color', [0,.8,0]);
   %p{jj,2} = plot(0:jplot-1,jr2_erf(jj,1:jplot), '-','linewidth',2, 'color', [0,0,.7]);
   
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
legend('RBC', 'JR');%, 'NK', 'NK-RBC', 'SW')

exportgraphics(f, 'figures_tables/tfp_compare_real.jpg','Resolution',300);


%% Additional statistics

return
%% Now make it a layered figure
for ss = 2:5
    for jj = 1:4
        p{jj,ss}.Visible = false;
    end
end

saveas(f, '../../../slides/WEAI_2021_alt/Figures/tfp_compare1.eps', 'epsc');

for ss = 2:5
   
    
    
    for jj = 1:4
        p{jj,ss-1}.Color = .75+min(.2*p{jj,ss-1}.Color,1); 
        p{jj,ss}.Visible = true;
    end
    
    saveas(f, ['../../../slides/WEAI_2021_alt/Figures/tfp_compare' num2str(ss) '.eps'], 'epsc');

end
