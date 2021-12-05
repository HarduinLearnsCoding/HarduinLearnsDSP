function [out] = myAbs32(vec)
%out=zeros(size(vec));
for k=1:length(vec) 
    out(k)=int32(abs(double(vec(k))));
end

end