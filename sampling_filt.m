function [k0] = sampling_filt(data,L)


for k=0:L-1
    E(k)=downsample(data(k+1:end-k),L);
end

[k0,emax]=max(E);

end
