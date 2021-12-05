function [out] = mySqr16(vec)
%out=zeros(size(vec));
for k=1:length(vec) 
    out(k)=MUL16_CX(CONJ32(vec(k)),vec(k));
end

end