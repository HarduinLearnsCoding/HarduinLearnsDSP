function [out]=sumall32(vec)
out=0;
for i=1:length(vec)
out=ADD32_CX(out,vec(i));
end
end