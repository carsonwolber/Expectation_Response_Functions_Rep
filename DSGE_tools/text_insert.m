function out = text_insert(fname_in,fname_out,text)

fid = fopen(fname_in, 'r');

fout = fopen(fname_out, 'w');




jj = 1;
while jj < 1000
    l = fgetl(fid);
    
    %Find all the bracketted pounds signs; will error if an unmatched pair
    idx = find(l=='#');
    while ~isempty(idx)
       l = [l(1:idx(1)-1), text{str2double(l(idx(1)+1:idx(2)-1))}, l(idx(2)+1:end)];
       idx = find(l=='#');
    end
    
    if l == -1
        break;
    end
    
    fprintf(fout,'%s\n', l);
    jj= jj+1;
end

fclose(fout);
fclose(fid);
