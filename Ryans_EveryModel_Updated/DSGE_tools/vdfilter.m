function [shares,variances] = vdfilter(A,B,C,lam,D)
% -------------------------------------------------------------------------
% Compute variance decomposition over frequencies lam = [lam1,lam2] when
% first differences have the state-space representation:
%       y(t) = a + A*x(t)   + D*u(t)
%       x(t) = b + B*x(t-1) + C*e(t)
% -------------------------------------------------------------------------

% Options
if nargin < 5; D = []; end;
if nargin < 4; lam = [2*pi/32,2*pi/6]; end; % default 6-32 periods

% Initialize
nx = size(B,1);
ne = size(C,2);
nu = size(D,2);
Vy = [];

% Compute filtered variances due to each shock
lam1 = lam(1);
lam2 = lam(2);
I    = eye(nx);


for i = 1:ne
    ind    = zeros(ne,1);
    ind(i) = 1;
    v      = C*ind;
%     
%     v = nsparse(v);
%     A = nsparse(A);
%     B = nsparse(B);
%     I = nsparse(I);
    
    %Old
%     fy     = @(w) 1/(2*pi)*(A/(I-B.*exp(-1i.*w))*v*(A/(I-B.*exp(1i.*w))*v).')./((1-exp(-1i.*w)).*(1-exp(1i.*w)));
%     Vy     = [Vy; diag(2.*real(integral(fy,lam1,lam2,'arrayvalued',true)))']    ;
%     

   %New
   %fy     = @(w) diag(1/(2*pi)*(A/(I-B.*exp(-1i.*w))*v*(A/(I-B.*exp(1i.*w))*v).')./((1-exp(-1i.*w)).*(1-exp(1i.*w))));
   
   fy = @(w) gomega(w,A,B,v,I);
   Vy     = [Vy; 2.*real(integral(fy,lam1,lam2,'arrayvalued',true))']    ;

end

% Include observation errors as shocks
if nargin > 4
    for i = 1:nu
        ind = zeros(nu,1);
        ind(i) = 1;
        v = D*ind;
        fy     = @(w) 1/(2*pi)*(v*v.')./((1-exp(-1i.*w)).*(1-exp(1i.*w)));
        Vy     = [Vy; diag(2.*real(integral(fy,lam1,lam2,'arrayvalued',true)))'];
    end
end
        
% Output objects
shares    = Vy./repmat(sum(Vy,1),ne+nu,1);
variances = sum(Vy,1); % filtered variances
end