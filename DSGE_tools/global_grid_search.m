% GLOBAL_GRID_SEARCH - Find the best value on an n^nparam grid. Depending
% on how fast your function is, can work for 2,3,or 4(??) variables.
%
% usage
%
% [param, min_val] = global_grid_search(f,lb,ub,n)

function [param, min_val] = global_grid_search(f,lb,ub,n)

dms = length(lb);
nel = n^dms;

%Make the grid
grids = cell(1,dms);
for jj = 1:dms
   grids{jj} = linspace(lb(jj),ub(jj),n); 
end
[grids{:}] = ndgrid(grids{:});
   
%Vectorize the grid
vals = zeros(dms,nel);
for jj = 1:dms
   vals(jj,:) = grids{jj}(:); 
end

%Compute objective at each point
out = zeros(1,nel);
parfor jj = 1:nel
    warning on
    out(jj) = f(vals(:,jj)');
end

%Find the best point
[min_val, v_idx] = min(out);
param = vals(:,v_idx)';