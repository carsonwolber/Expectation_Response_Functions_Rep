%**************************************************************
% SOLVE_LINEAR: Does some stuff.
%**************************************************************
addpath('DSGE_tools');

nper  = 20;
nerf  = 5000;

vnms = {'GDP', 'CON', 'PI', 'IR'};



%%
%**************************************************************
% NK ERF
%**************************************************************
clear *idx
[nk_param,nk_set] =nk_parameters;
load nk_obj

% Starndard parameters
nk_set.phipi =  1.5;       %taylor inflation
nk_set.phiy  =  0.5;       %taylor output level
nk_set.phidy =  0  ;       %taylor rule on output growth
nk_set.rhoi  =  0.5;       %taylor smoothing

[f, fx, fy, fxp, fyp, eta]=nk_prog(struct2array(nk_param),struct2array(nk_set));

yidx = [gap_idx:dy_idx];
xidx = [irl_idx:gapl_idx];
xkidx = ei_idx;
eqidx = 1:length([yidx,xidx]);


yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
std_erf = yerf([dy_idx,dy_idx,pi_idx,ir_idx],:);


%%
%**************************************************************
% BLL ERF
%**************************************************************
clear *idx
[bll_param,bll_set] =bll_parameters;
load bll_obj

% Starndard parameters
bll_set.phipi =  1.5;       %taylor inflation
bll_set.phiy  =  0.5;       %taylor output level
bll_set.phidy =  0  ;       %taylor rule on output growth
bll_set.rhoi  =  0.5;       %taylor smoothing

% %To get almost to RBC outcome...check why not more exact!
[f, fx, fy, fxp, fyp, eta]=bll_prog(struct2array(bll_param),struct2array(bll_set));

yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:al_idx];
xkidx = q_idx;
eqidx = 1:length([yidx,xidx]);
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
std_erf(:,:,3) = yerf([dy_idx,dc_idx,pi_idx,irt_idx],:);



%%
%**************************************************************
% NK ERF
%**************************************************************
clear *idx
[nk_param,nk_set] =nk_parameters;
load nk_obj

%BLL parameters
nk_set.phipi =  1.0137;       %taylor inflation
nk_set.phiy  = 0.0050;        %taylor output level
nk_set.phidy = 0;             %taylor rule on output growth
nk_set.rhoi  = 0.5583;        %taylor smoothing


[f, fx, fy, fxp, fyp, eta]=nk_prog(struct2array(nk_param),struct2array(nk_set));

yidx = [gap_idx:dy_idx];
xidx = [irl_idx:gapl_idx];
xkidx = ei_idx;
eqidx = 1:length([yidx,xidx]);

yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
bll_erf = yerf([dy_idx,dy_idx,pi_idx,ir_idx],:);


%%
%**************************************************************
% BLL ERF
%**************************************************************
clear *idx
[bll_param,bll_set] =bll_parameters;
load bll_obj

%BLL parameters
bll_set.phipi =  1.0137;       %taylor inflation
bll_set.phiy  = 0.0050;        %taylor output level
bll_set.phidy = 0;             %taylor rule on output growth
bll_set.rhoi  = 0.5583;        %taylor smoothing

% %To get almost to RBC outcome...check why not more exact!
[f, fx, fy, fxp, fyp, eta]=bll_prog(struct2array(bll_param),struct2array(bll_set));

yidx = [dy_idx:wpt_idx];
xidx = [kbar_idx:al_idx];
xkidx = q_idx;
eqidx = 1:length([yidx,xidx]);
yerf = gx_hx_erf(fy,fx,fyp,fxp,1,yidx,xidx,xkidx,eqidx,nerf);
bll_erf(:,:,3) = yerf([dy_idx,dc_idx,pi_idx,irt_idx],:);





%% Make Figure
f = figure;

p = cell(4,6);
for jj = 1:4
   s = subplot(2,2,jj);
    hold on 
   p{jj,1} = plot(0:jplot-1,std_erf(jj,1:jplot,1), '-b','linewidth',2);
   %p{jj,2} = plot(0:nerf,std_erf(jj,:,2), '-r','linewidth',2);
   p{jj,3} = plot(0:jplot-1,std_erf(jj,1:jplot,3), '-','linewidth',2, 'color',[.7,.5,1]);
   
   p{jj,4} = plot(0:jplot-1,bll_erf(jj,1:jplot,1), '-xb','linewidth',2);
   %p{jj,5} = plot(0:nerf,bll_erf(jj,:,2), '-xr','linewidth',2);
   p{jj,6} = plot(0:jplot-1,bll_erf(jj,1:jplot,3), '-x','linewidth',2, 'color',[.7,.5,1]);
   
   
   
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
legend('NK',  'Med. Scale', 'NK''',  'Med. Scale''')
legend('New Keyn.: stand. policy',  'Med. Scale: stand. policy', 'New Keyn.: accom. policy',  'Med. Scale: accom. policy', 'Location','NorthEast')


%% Now make it a layered figure
for ss = 4:6
    for jj = 1:4
        p{jj,ss}.Visible = 'off';
    end
end

saveas(f, 'figures_tables/fg_taylor_compare1.eps', 'epsc');

for ss = [1,3]
    for jj = 1:4
        p{jj,ss}.Color = .75+min(.2*p{jj,ss}.Color,1); 
        p{jj,ss+3}.Visible = 'on';
    end
end
saveas(f, 'figures_tables/fg_taylor_compare2.eps', 'epsc');

exportgraphics(f, 'figures_tables/fg_taylor_compare2.jpg','Resolution',300);
