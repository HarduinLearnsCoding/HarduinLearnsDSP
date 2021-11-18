function [x,y]= DFTdownsignal(xchannel,pilotsize)
xdown= downsample(xchannel,16);
%windowing the dft
for i=1:size(xdown)-pilotsize
    pilotpossibleloc=sum(abs(fft(xchannel[i:i+pilotsize])))
end
%checking max
pilotloc=max(pilotpossibleloc)



