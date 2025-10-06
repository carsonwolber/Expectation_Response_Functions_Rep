%**************************************************************
% SOLVE_LINEAR: Does some stuff.
%**************************************************************
addpath('../DSGE_tools');

nper =20;
load model_object
[param,set] = parameters;
[~,param] = model_ss(param,set);
paramv = struct2array(param);

%Get indexes to track
load v_idx
vlist     = {'GDP' 'PI',   'DA',   'DA+1|t',   'IR'};
vvidx     = [dy_idx,pi_idx, da_idx,dapt_idx,ir_idx];
prefs.ndec = 4;

%Define the neps shocks
shock = eye(neps); shcks = cell(1,neps);
for ss = 1:neps; shcks{ss} = shock(:,ss);end;
lnm = {'News-Perm', 'One-time Growth','Gov', 'Taylor', 'Noise', 'Alt Fund.', 'Alt. Noise' };
lstyle =  {'-b', '-r', '-g', '-m', '-k', '--r', '--k'};

%Cumulate growth-rates impulse responses to get level effects
cumsum_idx =[dy_idx,da_idx];

%**************************************************************
% FULL INFO MODEL USING BLL FUNDAMENTAL PROCESS
%**************************************************************
set.fullm = 1; set.bsinfo = 1;
[f, fx, fy, fxp, fyp, eta]=model_prog(paramv,struct2array(set));
[gx,hx]=gx_hx_alt(fy,fx,fyp,fxp);

disp(' ');
disp('FULL INFO VD:')
prefs.ndec = 4;
mom_tab(gx,hx,eta*eta',vvidx,vlist,prefs);disp(' ');
vdecom_table(gx,hx,eta,vvidx,vlist,prefs);disp(' ');
fevd_table(gx,hx,eta,[1,4,8,16,20],dy_idx,'DY',prefs);

%**************************************************************
% NOISE MODELS
%**************************************************************
forms = {'GenNoise','BS'};
set.bsinfo = 0;  % 0 = GenNoise, 1 = BS

set.fullm = 1; 
[f1, fx1, fy1, fxp1, fyp1, eta1]=model_prog(paramv,struct2array(set));
[gx1,hx1]=gx_hx_alt(fy1,fx1,fyp1,fxp1);

%Use version of model where PLM is parameter
set.fullm = 0; setv2 = struct2array(set);

%Exogenous impose gx1 PLM in full info forecasts
setv2(end-3*nx+1:end) = vec(gx1(e1t_idx:dapt_idx,:));

%Solve economy with gx1 PLM imposed
[f2, fx2, fy2, fxp2, fyp2, eta2,~, ss]=model_prog(paramv,setv2);
[gx2,hx2]=gx_hx_alt(fy2,fx2,fyp2,fxp2);

%Check
disp('CHECKING SANITY:');
disp(['Confirm PLM=ALM under full info- Residual = ' num2str(max(max(abs([gx2; hx2] - [gx1;hx1]))))])

%Pull out conjectured LOM
Gxmatn = gx2(e1t_idx:dapt_idx,:);

%For unknown states, move all dependence of beliefs from actual state
%to perceived state
if set.bsinfo
    plm_idx = (mu1t_idx:ggt_idx)-ny;
    alm_idx = (mu1_idx:gg_idx)-ny;
else
    plm_idx = (xxt_idx:evt_idx)-ny;
    alm_idx = (xx_idx:ev_idx)-ny;
end
Gxmatn(:,plm_idx) = Gxmatn(:,plm_idx)+Gxmatn(:,alm_idx);
Gxmatn(:,alm_idx) = 0;

%Re-solve the model
setv2(end-3*nx+1:end) = vec(Gxmatn);
[f2, fx2, fy2, fxp2, fyp2, eta2,~, ss]=model_prog(paramv,setv2);
[gx2,hx2]=gx_hx_alt(fy2,fx2,fyp2,fxp2);

disp(' ');
disp('NOISE REP:')
mom_tab(gx2,hx2,eta2*eta2',vvidx,vlist,prefs);disp(' ');
vdecom_table(gx2,hx2,eta2,vvidx,vlist,prefs);disp(' ');
fevd_table(gx2,hx2,eta2,[1,4,8,16,20],dy_idx,'DY',prefs);



clear mm simdata
eval(['save Solution_' forms{set.bsinfo+1} ' gx2 hx2 eta2 *_idx']);
