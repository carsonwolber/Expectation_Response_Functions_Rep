function out = nsparse(in, tol)

if nargin<2
    tol = 1e-10;
end

in(abs(in)<tol) = 0;

out = sparse(in);
