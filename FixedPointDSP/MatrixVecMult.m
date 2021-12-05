function [out]=MatrixVecMult(Mat,vec)
   
out=zeros(size(Mat,1),1);

for i=1:size(Mat,1)
for j=1:size(Mat,2)
    out(i)=ADD32_CX(out(i),MUL16_CX(Mat(i,j),vec(j)));
end
end

end

