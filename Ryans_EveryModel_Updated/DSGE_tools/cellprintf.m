function str = cellprintf(format,cell)
str = [];
for jj = 1:length(cell)
    str = [str, sprintf(format,cell{jj})];
end