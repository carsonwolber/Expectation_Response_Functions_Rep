function fmat = make_filter(n,bk_wghts,trend,s,S)
% computes (omega kron omega)*(LsS kron L0S)*F
% INPUTS:
%   n           = number of vars
%   bk_whgts    = BK-weights
%   trend       = location of trend-vars (empty for none)
%   s           = lag of cov that is to be computed (0 = contemp)
%   S           = lags of covariance vector that is multiplied

% infer k from BK_wghts
k  = (length(bk_wghts)-1)/2;

% make M1
M1 = zeros(1,n);
M1(trend)=1;
M1=diag(M1);

% computes sum of bk_weights(0:s-1).
beta_tilde = -fliplr(cumsum(bk_wghts));
%beta_tilde = -beta_tilde-bk_wghts; % flips order of x^k_-k, same as transposing Sig^k)

omega = kron(bk_wghts,eye(n)) + kron(beta_tilde,M1);

% quick-fix to implement first-difference filter
if k == 1
	beta_tilde = [0 1 0];
	omega = kron(bk_wghts,eye(n)-M1) + kron(beta_tilde,M1)
end


L0 = [speye(n*(2*k+1)), spzeros(n*(2*k+1),n*S)];
Ls = [spzeros(n*(2*k+1),n*s), speye(n*(2*k+1)), spzeros(n*(2*k+1),n*(S-s))];

OxO = kron(omega,omega);
LxL = kron(Ls,L0);
F   = construct_Fmat(n,2*k+S+1); 

fmat=OxO*LxL*F;