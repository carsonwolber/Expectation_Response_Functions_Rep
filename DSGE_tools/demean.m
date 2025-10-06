% DEMEAN - De-mean the columns of a matrix X
%
% usage
%
% out = demean(X)



function [out,mu] = demean(X)

t = size(X,1);
n = size(X,2);
mu = mean(X);

out = zeros(t,n);

for j = 1:n
    out(:,j) = X(:,j) - mu(j);
end

%out = X-ones(t,1)*mu;
