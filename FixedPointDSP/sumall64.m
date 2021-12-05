function [out]=sumall64(vec)
out=0;
for i=1:length(vec)
out=ADD64_CX(out,vec(i));
end
end