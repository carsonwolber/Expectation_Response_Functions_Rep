function out = pretty_print(V,ndec,rowtits)


sformat = ['%1.' num2str(ndec), 'f\t'];
for j = 1:size(V,1)
    str = [];
    if nargin>2
        str = sprintf('%s\t', rowtits{j});
    end
    disp([str, sprintf(sformat, V(j,:))]); 
end