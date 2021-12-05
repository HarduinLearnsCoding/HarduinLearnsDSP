function [Y] = myConvX(x,h)
m=length(x);
n=length(h);
X=[x;zeros(n,1)]; 
H=[h;zeros(m,1)]; 
for i=n:n+m-1
Y(i)=0;
for j=1:m
if(i-j+1>0)
Y(i)=ADD32_CX(Y(i),MUL16_CX(X(j),H(i-j+1)));
else
end
end
end

Y2=Y(n:end);

end