function out = rcumsum(x,dim)
disp('rcumsum \n');
if nargin == 1
    dim = 1;
end

if dim == 2
    x = x';
end

out = NaN(size(x));
for jj = 1:size(x,1)
    out(jj,:) = sum(x(jj:end,:));
end

if dim == 2
    out = out';
end

     

