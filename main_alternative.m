clear


% Sim paramaters
sim_set=1; % sim set 1 is SNR, sim set 2 is time delay and sim set 3 is freq delay
fs=16000;
f0=62.5;
t0=2.5/80;
p0=0.5; % phase error, from -0.5 to 0.5
SNR0=100;
FERlim=1000;
frame_number=1000;
% Fixed parameters
US=16; % Upsampling coeff
NumberofPilot=128;
NumberofKey=8;
NumberofDataSample=664;
NumberofAll=NumberofDataSample+NumberofKey+NumberofPilot;
BPSK_symbols=[-1;1];
%BPSK_symbols=[-1;1];
df=125;
dt=1/16;
th_f=2.5;

%LP Filter
n=[-45:45]; %sequence
pb=pi/16;
d=(pb)*sinc((pb).*n)-(-pb)*sinc((-pb).*n); %differences of 2 LP
hd = d.*hamming(length(d))'; %hamming for idealization

W128=ffttable(NumberofPilot,1,NumberofPilot,1,NumberofPilot);
W128trunc=ffttable(NumberofPilot,1,NumberofPilot,NumberofPilot/4/2-12+1,NumberofPilot/4/2+13);
% NOTE: As my calculations, the correct lowpass filter should be:
n=[-45:45]; 
pb=1/16;
d=(2*pb)*sinc((2*pb).*n); 
hd = d.*hamming(length(d))'; 
% However, the above shows more accurate performance 
hd=1; % All pass filter
hd = rcosdesign(0.35,6,16); % RCOS with 96 Taps

pi4tab=exp(1j*pi/4*[0:NumberofAll-1]');
switch sim_set
    case 0
        SNR_all=-3:0.5:15;
        it_num=length(SNR_all);
        f_all=f0*ones(1,it_num);
        t_all=t0*ones(1,it_num);
        p_all=p0*ones(1,it_num);
    case 1
        t_all=-2.5:0.5:2.5;
        %t_all=1/2*randi(10,1,30)-2.5;

        it_num=length(t_all);
        SNR_all=SNR0*ones(1,it_num);
        f_all=f0*ones(1,it_num);
        p_all=p0*ones(1,it_num);

    case 2
        f_all=125*randi(24,1,30)-1500;
        it_num=length(f_all);

        SNR_all=SNR0*ones(1,it_num);
        t_all=t0*ones(1,it_num);  
        p_all=p0*ones(1,it_num);
end


fprintf(1,'Simulation Status:      %3.1f%%',0);

for it=1:it_num
    
    p0_it=p_all(it);
    f0_it=f_all(it); 
    t0_it=t_all(it); 
    SNR_it=SNR_all(it);


    FrameError=0; % sum of Frame Errors 
    BitError=0;
    frame=0;
    errF=0;
    errT=0;
    errP=0;
    while (FrameError<FERlim && frame<frame_number)

        NP=10^(-SNR_it/10);
        
        data_ind=randi(2,NumberofDataSample,1);
        pilot_ind=2*ones(NumberofPilot,1);
        key_ind=ones(NumberofKey,1);
        data_BPSK=BPSK_symbols([pilot_ind;key_ind;data_ind]);
        data=data_BPSK.*exp(1j*pi/4*[0:NumberofAll-1]');
        
        % Upsampling + LP will come here
        data_us=upsample(data,US);
        data_lp=conv(data_us,hd,'same');
        %data_lp=lowpass(data_us,pi/16);

        data_dopp=data_lp.*(exp(2*pi*1j*f0_it/fs/US*[0:NumberofAll*US-1]'))*exp(2*pi*1j*p0_it);
        % Channel will come here
        data_tr=[zeros(US*40+t0_it/dt*US,1);data_dopp;zeros(US*40-t0_it/dt*US,1)];
        noise=sqrt(NP/2)*(randn(size(data_tr))+1j*randn(size(data_tr)));
        
        data_noisy=data_tr+noise;

        % Match Filter will come here
        data_lp2=conv(data_noisy,hd,'same');  
        
        % Sampling delay finder
        for k=0:US-1
            E1(k+1)=sum(mySqr(downsample(data_lp2(k+1:end-k),US)));
        end
        
        [emax,k0]=max(E1);
        k0=k0-1;
        data_ds=downsample(data_lp2(k0+1:end),US);
   
        
        for i=1:81 %CHANGED 10/29
            E2(i)=max(mySqr(W128*(data_ds(i:i+NumberofPilot-1))));
        end
        
        [~,nind]=max(E2);
        
        fft_plt=mySqr(W128trunc*(data_ds(nind:nind+NumberofPilot-1)));
%         After we have found the location of the pilot, we take the fft of
%         the pilot. We know that fft of the pilot signal has a peak in
%         2kHz when we transmitting it. After the channel block this peak
%         has shifted +/- 1500 Hz. Therefore, we expect a peak on the fft
%         of the pilot signal in the interval of 500Hz to 3500Hz. In fft it
%         is correspond to 5 to 29th indices. Multiplying with W128trunc
%         gives the these indices. Function 'mySqr' finds the square of the
%         magnitude of the vector.

         [m1,kind]=max(fft_plt);
         [Ks]=find(fft_plt>m1/th_f);
        [Ms]=fft_plt(Ks);
        
%         We are looking for peakes from the fft. Indice of the maximum
%         value gives us the where peak occurs. However, in order to get
%         more accurate solution, we also look the second peak. If there
%         exists some point larger than half of the peak value, we also
%         take this indice too.
        
        if length(Ks)>1 % 
            f_estimate=df*(sum(sqrt(Ms).*Ks)/(sum(sqrt(Ms)))-13); %CHANGED 10/29
        else
             f_estimate=df*(Ks-13);
        end
        
%          If there is exactly one peak, 'the indice number-13' multipy 125
%          Hz give us the frequency shift (Because when there is no shift, 
%          we have expect a peak in the 13th indice). On the other hand, if there are
%          more than one peak, we apply the linear interpolation. Sqrt
%          operation give us the magnitudes. Consider the magnitudes of peaks
%          are m1 and m2, and indices are k1 and k2, for first and second peak
%          respectively. Then linear interpolation can be written as k*=(m1k1+m2k2)/(m1+m2)

%         % Debug purpose
%         if not(nind==41)
%             1;
%         end
%         nind=41;
%         kind=1;
%         k0=0;

        % Debug Purpose
        
        t_estimate=dt/US*(nind*US+k0-3*US)-2.5;
        % 'Us*ks+nind' give us the in which indices our data has started.
        % dt is for indice to time conversion. If the data has start in the
        % first indice, it means there 2.5 time delay.

        % nindx=min(max(0,nind),880-NumberofAll); Not need anymore
        
        
        pilot_rec=exp(-2*pi*1j*f_estimate/fs.*[0:NumberofPilot-1]').*data_ds(nind:nind+NumberofPilot-1);
        pilot_dec=pilot_rec.*exp(-pi/4*1j*[0:NumberofPilot-1]');
        phase_dist=mean(pilot_dec(11:128))/BPSK_symbols(2); % because low pass filter, first indices of the pilot signal may lose some energy.
        % After we find, we first decode pilot data. Because we know
        % exactly what is our pilot data is, we can estimate phase shift
        % from that.
        
        p_estimate=atan2(imag(phase_dist),real(phase_dist))/pi/2;
        data_recieved=exp(-2*pi*1j*p_estimate)*exp(-2*pi*1j*f_estimate/fs.*[0:NumberofAll-1]').*data_ds(nind:nind+NumberofAll-1);
        % Now, we have estimated time shift, frequency shift and phase
        % delay, we can get the hat(Y) using all of these.
        
        % Symbol detection
        ind_recieved=sign(real(data_recieved.*exp(-1j*[0:NumberofAll-1]'*pi/4)))/2+1.5;
        % To demodulate the signal, we first convert BPSK pi/4 to BPSK,
        % then we apply the formula sing(x)/2+1.5 -> which gives 1 for
        % negative data, and 2 for positive data.
        
        data_ind_received=ind_recieved(NumberofPilot+NumberofKey+1:end);
        % We can get indices correspond to the data from the demodulated signal in order to count error. 
        
        BitErrorCount=sum(nnz(data_ind-data_ind_received));
        % Compare received indices with transmitted indices, to find BER.
        
        FrameErrorCount=BitErrorCount>0;
        % Look for frame error.
        
        BitError=BitError+BitErrorCount;
        FrameError=FrameError+FrameErrorCount;
        frame=frame+1;
        errF(frame)=f_estimate-f0_it;
        errT(frame)=t_estimate-t0_it;
        errP(frame)=p_estimate-p0_it;

   
        
    end
    BER(it)=BitError/(NumberofDataSample*frame);
    FER(it)=FrameError/frame;
    mseF(it)=sqrt(mean(errF.^2));
    mseT(it)=sqrt(mean(errT.^2));
    mseP(it)=sqrt(mean(errP.^2));
    frames(it)=frame;
    fprintf(1,'\b\b\b\b\b%3.1f%%',100*it/(it_num));
end
        dmin=mySqr(BPSK_symbols(1)-BPSK_symbols(2));
        e=10.^(SNR_all/10);        
        theoric = qfunc(sqrt((dmin^2)*(e/2)));

% Figures
switch sim_set
    case 0
        dmin=abs(BPSK_symbols(1)-BPSK_symbols(2));
        e=10.^(SNR_all/10);        
        theoric = qfunc(sqrt((dmin^2)*(e/2)));
        th=qfunc(sqrt(2*e));

    figure;
    semilogy(SNR_all,BER,'-x')
    hold on;
    %semilogy(SNR_all,th,'+')
    semilogy(SNR_all,theoric,'--o')
    grid on;
    xlim([min(SNR_all),max(SNR_all)])
    %ylim([min(BER(BER>0)),1])
    ylim([1e-4,1])
    xlabel('SNR in dB')
    ylabel('BER')
    legend('Simulation Result','Ideal BPSK BER Curve')
    
    
     figure;
    semilogy(SNR_all,FER)
    grid on;
    xlim([min(SNR_all),max(SNR_all)])
    %ylim([min(BER(BER>0)),1])
    %ylim([1e-4,1])
    xlabel('SNR in dB')
    ylabel('FER')
    
        figure;
    plot(SNR_all,mseF)
    grid on;
        xlabel('SNR in dB')
    ylabel('STD of Frequency Estimation')
        xlim([min(SNR_all),max(SNR_all)])

          figure;
    plot(SNR_all,mseT)
    grid on;
        xlabel('SNR in dB')
    ylabel('STD of Time Estimation')
        xlim([min(SNR_all),max(SNR_all)])

    case 1
        
        
        
  figure;
    plot(t_all,'-x')
    hold on;
    plot(mseT+t_all,'o')
    grid on;
    xlabel('Number of trial')
    ylabel('Time (s)')
    %xlabel('Frequency shift (Hz)')
    %ylabel('STD of Frequency Estimation')
    legend('Actual Time Shift','Time Shift Estimation')
        
    
   
      figure;
    plot(t_all,mseT)
    grid on;
    xlabel('Time Delay (ms)')
    ylabel('STD of Time Estimation')
    
            figure;
    plot(f_all,mseF)
    grid on;
    xlabel('Frequency shift (Hz)')
    ylabel('STD of Frequency Estimation')
    
        figure;
    semilogy(t_all,BER)
    hold on;
    semilogy(t_all,theoric,'o')
    grid on;
    xlim([min(t_all),max(t_all)])
    ylabel('BER')
    xlabel('SNR in dB')
    
         figure;
    semilogy(t_all,FER)
    grid on;
    xlim([min(t_all),max(t_all)])
    %ylim([min(BER(BER>0)),1])
    %ylim([1e-4,1])
    xlabel('SNR in dB')
    ylabel('FER')
    case 2
        figure;
    plot(f_all,mseF)
    grid on;
    xlabel('Frequency shift (Hz)')
    ylabel('STD of Frequency Estimation')
        xlim([min(f_all),max(f_all)])

  figure;
    plot(f_all,'-x')
    hold on;
    plot(mseF+f_all,'o')
    grid on;
    xlabel('Number of trial')
    ylabel('Frequency (Hz)')
    %xlabel('Frequency shift (Hz)')
    %ylabel('STD of Frequency Estimation')
    legend('Actual Freq Shift','Freq Shift Estimation')
        
    
          figure;
    plot(f_all,mseT)
    grid on;
    xlabel('Frequency shift (Hz)')
    ylabel('STD of Time Estimation')
    
            figure;
    semilogy(f_all,BER)
    hold on;
    semilogy(f_all,theoric,'o')
    grid on;
    xlim([min(f_all),max(f_all)])
    ylabel('BER')
        xlabel('SNR in dB')
        
        
         figure;
    semilogy(f_all,FER)
    grid on;
    %ylim([min(BER(BER>0)),1])
    %ylim([1e-4,1])
    xlabel('SNR in dB')
    ylabel('FER')
end


