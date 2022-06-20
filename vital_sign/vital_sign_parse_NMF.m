clear;
[retVal] = readDCA1000_1('C:\Users\Administrator\Desktop\�״�����\0.3-1wr.bin');
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

%NMF
[w,h]=nmf(RX1data,5,200);
nmf_result=w*h;

%�Ӻ�����
range_win = hamming(numADCSamples); %���ɺ�����
for k=1:1:frame
    din_win(:,k)=nmf_result(:,2*k-1).*range_win; %���ź���Range-fft
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

%��ȡ�ź�ʵ�����鲿
for k=1:frame
    data_real(:,k)=real(data(:,k));
    data_imag(:,k)=imag(data(:,k));
end

%�����ź���λ
for k=1:frame
    signal_phase(:,k)=atan(data_imag(:,k)/data_real(:,k));
end

%��λչ��
for k=2:frame
    diff=signal_phase(:,k)-signal_phase(:,k-1);
    if diff>pi/2
        signal_phase(:,(k:end))=signal_phase(:,(k:end))-pi;
    elseif diff<-pi/2
        signal_phase(:,(k:end))=signal_phase(:,(k:end))+pi;
    end
end

%������λ��
for k=1:frame-1
    delta_phase(:,k)=signal_phase(:,k+1)-signal_phase(:,k);
end

%�Ӳ�����ȥ����������
thresh=0.8;
for k=1:frame-3
    phaseUsedComputation(:,k)=filter_RemoveImpulseNoise(delta_phase(:,k),delta_phase(:,k+1),delta_phase(:,k+2),thresh);
end
index=1:1:frame-3;
index=index*Tf;

%ԭʼ�źŴ�ͨ�˲�
filter_delta_phase=filter(bpf_vitalsign,phaseUsedComputation);

%����λת��Ϊλ��
% delta_d=(filter_delta_phase*lambda)/4/pi;

%��ʵ��ǻλ���ź�
% for i=1:frame-3
%    vital_sign(i)= (filter_delta_phase(1)+filter_delta_phase(i))*i/2;
% end

vital_sign=filter_delta_phase;

%ԭʼ�ź�ʱ��ͼ
figure(1);
plot(index,vital_sign);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('�ķ��ź�','FontWeight','bold');

%��ԭʼ�ź���fft
vital_sign_fft=fft(vital_sign,N);

%˫�ߴ��ź�תΪ���ߴ�
freq=(0:1:N/2)/Tf/N;  %���������źŲ�����Ϊ֡��
P2 = abs(vital_sign_fft/(N-1));
P1 = P2(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1(2:end-1) = 2*P1(2:end-1);

%ԭʼ�ź�Ƶ��ͼ
figure(2);
plot(freq,P1);
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('�ķ��ź�Ƶ��ͼ','FontWeight','bold');

%�����źŴ�ͨ�˲�
filter_delta_phase_breathe=filter(bpf_breathe,phaseUsedComputation);

%�ź�תΪ��ʵ�����ź�
% for i=1:frame-3
%    breathe(i)= (filter_delta_phase_breathe(1)+filter_delta_phase_breathe(i))*i/2;
% end

breathe=filter_delta_phase_breathe;

%�����ź�ʱ��ͼ
figure(3);
plot(index,breathe);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('�����ź�','FontWeight','bold');

%�Ժ����ź���fft
breathe_fft=fft(breathe,N);

%˫�ߴ��ź�תΪ���ߴ�
P2_breathe = abs(breathe_fft/(N-1));
P1_breathe = P2_breathe(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1_breathe(2:end-1) = 2*P1_breathe(2:end-1);

%�����ź�Ƶ��ͼ
figure(4);
plot(freq,P1_breathe);
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('�����ź�Ƶ��ͼ','FontWeight','bold');


%�����źŴ�ͨ�˲�
filter_delta_phase_heart=filter(bpf_heart,phaseUsedComputation);

%��ʵ����λ���ź�
% for i=1:frame-3
%    heart(i)= (filter_delta_phase_heart(1)+filter_delta_phase_heart(i))*i/2;
% end

heart=filter_delta_phase_heart;

%�����ź�ʱ��ͼ
figure(5);
plot(index,heart);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold ');
title('�����ź�','FontWeight','bold');

%�������ź���fft
heart_fft=fft(heart,N);

%˫�ߴ��ź�תΪ���ߴ�
P2_heart = abs(heart_fft/(N-1));
P1_heart = P2_heart(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1_heart(2:end-1) = 2*P1_heart(2:end-1);

%����г�����
[heart_peaks,heart_peaksnum]=findpeaks(P1_heart,0.9,2,N,Tf);%0.9Hz-2Hz
[heart_harmonic_peaks,heart_harmonic_peaksnum]=findpeaks(P1_heart,1.8,4,N,Tf);%1.8Hz-4Hz
heart_peaks=heart_peaks/N/Tf;
heart_harmonic_peaks=heart_harmonic_peaks/N/Tf;
[heart_peaks_row,heart_peaks_column]=size(heart_peaks);
[heart_harmonic_peaks_row,heart_harmonic_peaks_column]=size(heart_harmonic_peaks);
for i=1:heart_peaks_column
    if max(P1_heart)-P1_heart(round(heart_peaks(i)*N*Tf)+1)<0.3 %����������ֵ֮����һ����Χ�ڣ��ٽ��ҵ�г�����źŽ�������
        for j=1:heart_harmonic_peaks_column
            if heart_harmonic_peaks(j)/heart_peaks(i)==2
                P1_heart(round(heart_peaks(i)*N*Tf)+1)=2*P1_heart(round(heart_peaks(i)*N*Tf)+1);
            end
        end
    end
end

%�����ź�Ƶ��ͼ
figure(6);
plot(freq,P1_heart);
xlim([0,4]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('�����ź�Ƶ��ͼ','FontWeight','bold');