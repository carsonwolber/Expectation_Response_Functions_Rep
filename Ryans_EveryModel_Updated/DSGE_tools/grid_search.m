%grid_search(f,start,lb,ub,n)

function [param, min_val] = grid_search(f,start,lb,ub,n)

if length(start) ~= length(lb)
    error('length(start) != length(lb)');
end

%Make grid
grid = zeros(length(lb),n+1);
for j = 1:length(lb)
    grid(j,:) = [linspace(lb(j)+.01*abs(lb(j)),ub(j)-.01*abs(ub(j)),n),start(j)];
end


%Initialize values
idx = (n+1) + zeros(length(lb),1);
min_val = f(start);
min_new = min_val;
param = start;


i = 1;
return_id = false;
disp(['Intial fval: ', num2str(min_val)]);
f_new = zeros(1,n+1);
while i <8 && ~return_id
    for k = 1:length(lb)
        parfor j = 1:n+1
            param_tmp = param;
            param_tmp(k) = grid(k,j);
            f_new(j) = f(param_tmp);
        end
        [min_tmp,idx_tmp] = min(f_new);
        
        if min_tmp<min_new
            param(k) = grid(k,idx_tmp);
            min_new = min_tmp;
        end
    end
    
    %Update Things for nex loop
    disp(['Improvement: ', num2str(min_val-min_new)]);
    if min_new==min_val
        return_id = true;
    end
    min_val = min_new;
    disp(['fval: ', num2str(min_val)])
    disp(' ');
    i = i+1;

end




