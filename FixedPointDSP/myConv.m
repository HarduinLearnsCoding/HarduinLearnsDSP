function [out] = myConv(v1,v2)
L1=length(v1);
L2=length(v2);
out=zeros(size(v1));
for k=1:L1
    for s=1:max(L1,L2)
        if ((k-s+1)>0 && (k-s+1)<L2+1)
        out(k)=out(k)+v1(s)*v2(k-s+1);
        end
    end
end
end