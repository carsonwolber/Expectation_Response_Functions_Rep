%Function out = selector(idx,N)

function out = selector(idx,N)

out = zeros(length(idx),N);

for jj = 1:length(idx)
    
    out(jj,idx(jj)) = 1;
end