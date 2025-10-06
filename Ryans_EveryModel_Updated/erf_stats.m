function [SFW,MHE,LPF32,LPF100] = erf_stats(erf_in,GAM,ZZ,TQA,TQB)


%Share future weights
%SFW = sum(abs(erf_in(:,2:end,:).^2),2)./sum(abs(erf_in.^2),2);


SFW = diag(erf_in(:,2:end,:)*erf_in(:,2:end,:)')./diag(erf_in(:,1:end,:)*erf_in(:,1:end,:)');


%Mean horizon of expectation
%MHE = sum(abs(erf_in.^2).*(0:size(erf_in,2)-1),2)./sum(abs(erf_in.^2),2)

MHE = diag((0:(size(erf_in,2)-1)).*(erf_in)*erf_in')./diag(erf_in(:,1:end,:)*erf_in(:,1:end,:)');

%FD measures
v_all   = diag(erf_filter(GAM,ZZ,TQA,TQB,[0*2*pi/10000000,2*pi/2]));
v_lp32  = diag(erf_filter(GAM,ZZ,TQA,TQB,[0*2*pi/10000000,2*pi/32]));
v_lp100 = diag(erf_filter(GAM,ZZ,TQA,TQB,[0*2*pi/10000000,2*pi/100]));


LPF32  = v_lp32./v_all;
LPF100 = v_lp100./v_all;


