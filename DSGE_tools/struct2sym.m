function PARAM = struct2sym(param)



param_list = fieldnames(param);
maxp = inf;
PARAM = sym(zeros(1,length(param_list)));
for j = 1:min(length(param_list),maxp)
    eval([param_list{j} ' =sym(param_list{j});']);
    assignin('caller', param_list{j}, eval(param_list{j}));
    eval(['PARAM(j) = ',param_list{j} ,';']);
    if j==maxp
        disp('Max sym parameters reached');
    end
end
