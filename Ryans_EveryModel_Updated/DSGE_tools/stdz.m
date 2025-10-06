%Standize data to have mean 0 and variance xx by cols

function [out, weights, mu] = stdz(X, varargin)

t = size(X,1);
n = size(X,2);
mu = mean(X);

out = X-ones(t,1)*mu;

weights = sqrt(sum((out.^2)))/sqrt((t-1));

if ~isempty(varargin) && strcmp(varargin{1},'std*(t-1)')
    %weights = std(out)*sqrt(t-1);
    weights = weights*sqrt(t-1);
    
else
    %weights = std(out);
    %weights2 = sqrt(sum((out.^2)))/sqrt((t-1));
end

out = out./(ones(t,1)*weights);
