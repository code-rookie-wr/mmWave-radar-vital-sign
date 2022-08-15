% ���ļ����ڶ��ڳ�ʱ����������Լ����ʵļ�⣬��demo����20s���ݣ���������ο�
clear;
[retVal] = readDCA1000_1('./demo.bin');
global numChirps;
global numADCSamples;%��������
RX1data = reshape(retVal(1,:),numADCSamples,numChirps);   %RX1����
RX2data = reshape(retVal(2,:),numADCSamples,numChirps);   %RX2
RX3data = reshape(retVal(3,:),numADCSamples,numChirps);   %RX3
RX4data = reshape(retVal(4,:),numADCSamples,numChirps);   %RX4

c=3.0e8;  
slope=60e12; %��Ƶб��
Tc=50e-6;      %chirp����
B=slope*Tc;    %��Ƶ����
Fs=4e6;        %������
f0=60.36e9;    %��ʼƵ��
lambda=c/f0;   %�״��źŲ���
d=lambda/2;    %�������м��
frame=400;     %֡��
Tf=0.05;       %֡����
N=1024;        %FFT����

range_win = hamming(numADCSamples); %���ɺ�����
for k=1:1:frame
    din_win(:,k)=RX1data(:,2*k-1).*range_win; %���ź���Range-fft
    datafft(:,k)=fft(din_win(:,k));
end

%�ҵ�Range-bin��ֵ
rangeBinStartIndex=3;%�ֱ���0.1m
rangeBinEndIndex=10;
for k=1:1:frame
    for j=rangeBinStartIndex:1:numADCSamples 
        if(abs(datafft(j,k))==max(abs(datafft((rangeBinStartIndex:rangeBinEndIndex),k)))) 
            data(:,k)=datafft(j,k);
        end
    end
end

%�����ݼӸ�300֡���ڣ�ÿ����20֡���л���
frame_process=300;
slide_width=20;
for i=0:((frame-frame_process)/slide_width)
    data_process(i+1,:)=data(:,(1+slide_width*i):(frame_process+slide_width*i));

    
%��ȡ�ź�ʵ�����鲿
for k=1:frame_process
    data_real(i+1,k)=real(data_process(i+1,k));
    data_imag(i+1,k)=imag(data_process(i+1,k));
end

%�����ź���λ
for k=1:frame_process
    signal_phase(i+1,k)=atan(data_imag(i+1,k)/data_real(i+1,k));
end

%��λչ��
for k=2:frame_process
    diff=signal_phase(i+1,k)-signal_phase(i+1,k-1);
    if diff>pi/2
        signal_phase(i+1,(k:end))=signal_phase(i+1,(k:end))-pi;
    elseif diff<-pi/2
        signal_phase(i+1,(k:end))=signal_phase(i+1,(k:end))+pi;
    end
end

%������λ��
for k=1:frame_process-1
    delta_phase(i+1,k)=signal_phase(i+1,k+1)-signal_phase(i+1,k);
end

%�Ӳ�����ȥ����������
thresh=1;
for k=1:frame_process-3
    phaseUsedComputation(i+1,k)=filter_RemoveImpulseNoise(delta_phase(i+1,k),delta_phase(i+1,k+1),delta_phase(i+1,k+2),thresh);
end

%�����źŴ�ͨ�˲�
filter_delta_phase_breathe(i+1,:)=filter(bpf_breathe,phaseUsedComputation(i+1,:));

%�ź�תΪ��ʵ�����ź�
% for k=1:frame_process-3
%    breathe(i+1,k)= (filter_delta_phase_breathe(1)+filter_delta_phase_breathe(i+1,k))*k/2;
% end

breathe=filter_delta_phase_breathe;

%�Ժ����ź���fft
breathe_fft(i+1,:)=fft(breathe(i+1,:),N);

%˫�ߴ��ź�תΪ���ߴ�
P2_breathe = abs(breathe_fft(i+1,:)/(N-1));
P1_breathe = P2_breathe(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1_breathe(2:end-1) = 2*P1_breathe(2:end-1);
ssb_breathe(i+1,:)=P1_breathe;

%�ҵ�����Ƶ���ֵ
breathe_max(i+1)=findpeaksmax(ssb_breathe(i+1,:),0.1,0.6,N,Tf)/N/Tf*60;

%�����źŴ�ͨ�˲�
filter_delta_phase_heart(i+1,:)=filter(bpf_heart,phaseUsedComputation(i+1,:));

%��ʵ����λ���ź�
% for k=1:frame_process-3
%    heart(i+1,k)= (filter_delta_phase_heart(1)+filter_delta_phase_heart(i+1,k))*k/2;
% end

heart=filter_delta_phase_heart;

%�������ź���fft
heart_fft((i+1),:)=fft(heart((i+1),:),N);

%˫�ߴ��ź�תΪ���ߴ�
P2_heart = abs(heart_fft((i+1),:)/(N-1));
P1_heart = P2_heart(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1_heart(2:end-1) = 2*P1_heart(2:end-1);

%����г�����
[heart_peaks,heart_peaksnum]=findpeaks(P1_heart,0.9,2,N,Tf);%0.9Hz-2Hz
[heart_harmonic_peaks,heart_harmonic_peaksnum]=findpeaks(P1_heart,1.8,4,N,Tf);%1.8Hz-4Hz
heart_peaks=heart_peaks/N/Tf;
heart_harmonic_peaks=heart_harmonic_peaks/N/Tf;
[heart_peaks_row,heart_peaks_column]=size(heart_peaks);
[heart_harmonic_peaks_row,heart_harmonic_peaks_column]=size(heart_harmonic_peaks);
for k=1:heart_peaks_column
    if max(P1_heart)-P1_heart(round(heart_peaks(k)*N*Tf)+1)<0.3 %����������ֵ֮����һ����Χ�ڣ��ٽ��ҵ�г�����źŽ�������
        for j=1:heart_harmonic_peaks_column
            if heart_harmonic_peaks(j)/heart_peaks(k)==2
                P1_heart(round(heart_peaks(k)*N*Tf)+1)=2*P1_heart(round(heart_peaks(k)*N*Tf)+1);
            end
        end
    end
end
ssb_heart(i+1,:)=P1_heart;

%�ҵ�����Ƶ���ֵ
heart_max(i+1)=findpeaksmax(ssb_heart(i+1,:),0.9,2,N,Tf)/N/Tf*60;
end

%�趨ʱ��ָ��
index=(frame_process/20):(slide_width/20):(frame/20);

%���ƺ�����̬�仯ͼ
figure(1);
plot(index,breathe_max);
xlabel('time(s)');
ylabel('times/min');
title('����������̬�仯ͼ');

%����������̬�仯ͼ
figure(2);
plot(index,heart_max);
xlabel('time(s)');
ylabel('times/min');
title('����������̬�仯ͼ');