function [out] = myAbs16(vec)
%out=zeros(size(vec));
for k=1:length(vec) 
    out(k)=int16(abs(double(vec(k))));
end

end