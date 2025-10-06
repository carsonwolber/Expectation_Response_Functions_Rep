function out = table_insert(fin,fout,data,sformats)

sout = cell(size(data));
for rr = 1:size(data,1)
    for cc = 1:size(data,2)
        if ~isnan(data(rr,cc))
            sout{rr,cc} = sprintf(sformats{min(length(sformats),rr)},data(rr,cc));
        else
            sout{rr,cc} = '-';
        end
    end
end
sout = transpose(sout);

text_insert(fin,fout,sout(:));