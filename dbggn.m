clear


% %LP Filter
% n=[-45:45]; %sequence
% pb=pi/16;
% d=(pb)*sinc((pb).*n)-(-pb)*sinc((-pb).*n); %differences of 2 LP
% hd = d.*hamming(length(d))'/1; %hamming for idealization
% 
% 
% n=[-50:50];
% 
% bp=[697,770,852,941,1209,1336,1477,1633]./(fs/2);
% pb=30/fs;
% for i=1:8
% d(i,:)=(bp(i)+pb)*sinc((bp(i)+pb).*n)-(bp(i)-pb)*sinc((bp(i)-pb).*n);
% d(i,:) = d(i,:).*hamming(101)';
% end



% Sim paramaters
sim_set=0; % sim set 1 is SNR, sim set 2 is time delay and sim set 3 is freq delay
fs=16000;
f0=0;
t0=0;
SNR0=100;
FERlim=100;
frame_number=500;
US=16; % Upsampling coeff
% Fixed parameters
NumberofPilot=128;
NumberofKey=8;
NumberofDataSample=664;
NumberofAll=NumberofDataSample+NumberofKey+NumberofPilot;
BPSK_symbols=[-1-1j;1+1j]/sqrt(2);
%BPSK_symbols=[-1;1];
df=125;
dt=1/16;


n=[-45:45];
bp=1000/fs/2;

d=(2*bp)*sinc((2*bp).*n);
hd_sinc1 = d.*hamming(length(d))';


bp=pi/32;

d2=(bp)*sinc((bp).*n);
hd_sinc2 = d2.*hamming(length(d2))';
% fvtool(hd)
%hd=1;
% % 
N = 90;
Fs = 16000;
Fp = 1000;
Ap = 0.01;
Ast = 80;
Rp = (10^(Ap/20) - 1)/(10^(Ap/20) + 1);
Rst = 10^(-Ast/20);
hd_firce = firceqrip(N,Fp/Fs,[Rp Rst],'passedge');

% n=[-45:45]; 
% pb=1/16;
% d=(2*pb)*sinc((2*pb).*n); 
% hd = d.*hamming(length(d))'; 
%  fvtool(NUM)
% 
% fvtool(d,'Fs',fs)
hd_no=[1];

hd=hd_sinc2;
switch sim_set
    case 0
        SNR_all=-3:5;
        it_num=length(SNR_all);
        f_all=f0*ones(1,it_num);
        t_all=t0*ones(1,it_num);
    case 1
        t_all=-2.5:0.5:2.5;
        it_num=length(t_all);
        SNR_all=SNR0*ones(1,it_num);
        f_all=f0*ones(1,it_num);
    case 2
        f_all=0:125:2000;
        it_num=length(f_all);
        SNR_all=SNR0*ones(1,it_num);
        t_all=t0*ones(1,it_num);    
end


fprintf(1,'Simulation Status:      %3.1f%%',0);

for it=1:it_num
    
    f0_it=f_all(it); 
    t0_it=t_all(it); 
    SNR_it=SNR_all(it);


    FrameError=0; % sum of Frame Errors 
    BitError=0;
    frame=0;
    errF=0;
    errT=0;
    while (FrameError<FERlim && frame<frame_number)

        NP=10^(-SNR_it/10);
        
        data_ind=randi(2,NumberofDataSample,1);
        pilot_ind=2*ones(NumberofPilot,1);
        key_ind=ones(NumberofKey,1);
        data=BPSK_symbols([pilot_ind;key_ind;data_ind]);
        
        % Upsampling + LP will come here
        data_us=upsample(data,US);
        data_lp=conv(data_us,hd,'same');
        %data_lp=lowpass(data_us,pi/16);

        data_dopp=data_lp.*(exp(2*pi*1j*f0_it/fs*[0:NumberofAll*US-1]'));
        % Channel will come here
        data_tr=[zeros(US*40+t0_it/dt*US,1);data_dopp;zeros(US*40-t0_it/dt*US,1)];
        noise=sqrt(NP/2)*(randn(size(data_tr))+1j*randn(size(data_tr)));
        
        data_noisy=data_tr+noise;

        % Match Filter will come here
          data_lp2=conv(data_noisy,hd,'same');  
          %data_lp2=data_noisy;
          %data_lp2=lowpass(data_noisy,pi/16);
        % Sampling delay finder
        %[k0]=sampling_filt(data_lp,US);
        for k=0:US-1
            E1(k+1)=sum(abs(downsample(data_lp2(k+1:end-k),US)).^2);
        end
        
        [emax,k0]=max(E1);
        k0=k0-1;
        data_ds=downsample(data_lp2-k0,US);
   
        for i=1:length(data_ds)-NumberofPilot
            E2(i)=max(abs(fft(data_ds(i:i+NumberofPilot-1))));
%             F(:,i)=abs(fft(data_ds(i:i+NumberofPilot-1)));
%             %figure(2)
%             %plot(F)
        end
        
        [~,nind]=max(E2);
        [~,kind]=max(abs(fft(data_ds(nind:nind+NumberofPilot-1))));
        
        % Debug purpose
%         if not(nind==41)
%             1;
%         end
%         nind=41;
%         kind=1;
%         k0=0;
        
        
        f_estimate=df*(kind-1);
        t_estimate=dt*(nind+k0/US-1)-2.5;
        
        nindx=min(max(0,nind),880-NumberofAll);
        
        data_recieved=exp(-2*pi*f_estimate/fs*1j.*[0:US:NumberofAll*US-1])'.*data_ds(nindx:nindx+NumberofAll-1);
        % Symbol detection
        ind_recieved=sign(real(data_recieved*exp(-1j*pi/4)))/2+1.5;
        data_ind_received=ind_recieved(NumberofPilot+NumberofKey+1:end);
        BitErrorCount=sum(nnz(data_ind-data_ind_received));
        FrameErrorCount=BitErrorCount>0;
        BitError=BitError+BitErrorCount;
        FrameError=FrameError+FrameErrorCount;
        frame=frame+1;
        errF(frame)=f_estimate-f0_it;
        errT(frame)=t_estimate-t0_it;
        
        if (BitErrorCount/NumberofDataSample>0.4)
        1;
        end
        
    end
    BER(it)=BitError/(NumberofDataSample*frame);
    FER(it)=FrameError/frame;
    mseF(it)=mean(errF.^2);
    mseT(it)=mean(errT.^2);

    fprintf(1,'\b\b\b\b\b%3.1f%%',100*it/(it_num));
end


% Figures
switch sim_set
    case 0
        dmin=abs(BPSK_symbols(1)-BPSK_symbols(2));
        e=10.^(SNR_all/10);        
        theoric = qfunc(sqrt((dmin^2)*(e/2)));
        th=qfunc(sqrt(2*e));

    figure;
    semilogy(SNR_all,BER)
    hold on;
    semilogy(SNR_all,th,'+')
    semilogy(SNR_all,th,'o')
    grid on;
    
%     figure;
%     semilogy(SNR_all,FER)
    case 1
    mse_of_timedelay=(mseT).^2;
      figure;
    plot(t_all,mseT)
    grid on;
    xlabel('Time Delay (ms)')
    ylabel('MSE of Time Estimation')
    case 2
        figure;
    plot(t_all,mseF)
    grid on;
    xlabel('Frequency shift (Hz)')
    ylabel('MSE of Frequency Estimation')
end


