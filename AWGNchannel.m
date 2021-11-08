%Creating the channel


function [xafterchannel] = AWGNchannel(xmod,EbNodb,SNRdb)

% powxmod= sum(abs(xmod).^2)/length(xmod);

EbNo=10^(EbNodb/10);
SNRdb=10^(SNRdb/10);
N0=Eb/SNR
wgn=(randn(1,length(xmod))+1i*randn(1,length(xmod)))*sqrt(No/2);
xafterchannel= xmod + wgn;


f0= randi([-1500 1500],columns,1); 
k0= randi([-64 64],1, columns);
r=[];
noisechannel=[];
r(timenew)= exp(1i*2*pi*f0*timenew)*(timenew-k0) + noisechannel(timenew);



