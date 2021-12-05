function [out]=sumall16(vec)
out=0;
for i=1:length(vec)
out=ADD16_CX(out,vec(i));
end
end