function [dftblock]=ffttable(size,column_start,column_end,row_start,row_end)
dftblock=zeros(size,size);
for i=1:size
    for j=1:size
        dftblock(i,j)=exp(-(1i*2*pi*(i-1)*(j-1))/size);
    end
end
dftblock=dftblock(row_start:row_end,column_start:column_end);
end
