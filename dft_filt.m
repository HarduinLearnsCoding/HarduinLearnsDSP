function [n0,f0] = dft_filt(data,psize,US)

Len=length(data);

for i=1:Len-psize
    E(i)=sum(abs(fft(data(i:i+psize-1))));
end

[nind,~]=max(E);
[kind,~]=max(abs(fft(data(nind:nind+psize-1))));



end

