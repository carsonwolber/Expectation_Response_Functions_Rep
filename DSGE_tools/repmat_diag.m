function out = repmat_diag(X,n)

out = kron(speye(n),X);