clear;
[retVal] = readDCA1000_1('./demo.bin');
global numChirps;
global numADCSamples;%��������
RX1data = reshape(retVal(1,:),numADCSamples,numChirps);   %RX1����
RX2data = reshape(retVal(2,:),numADCSamples,numChirps);   %RX2
RX3data = reshape(retVal(3,:),numADCSamples,numChirps);   %RX3
RX4data = reshape(retVal(4,:),numADCSamples,numChirps);   %RX4

c=3.0e8;  
slope=59.997e12; %��Ƶб��
Tc=50e-6;      %chirp����
B=slope*Tc;    %��Ƶ����
Fs=4e6;        %������
f0=60.36e9;    %��ʼƵ��
lambda=c/f0;   %�״��źŲ���
d=lambda/2;    %�������м��
frame=400;     %֡��
Tf=0.05;       %֡����
N=1024;        %FFT����

%�Ӻ�����
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

%ԭʼ�źŴ�ͨ�˲�
filter_delta_phase=filter(bpf_vitalsign,phaseUsedComputation);

vital_sign=filter_delta_phase;

%������ʱ������
for k=4:frame-6
   vital_sign_diff1(k-3)=(5*(vital_sign(k+1)-vital_sign(k-1))+4*(vital_sign(k+2)-vital_sign(k-2))+(vital_sign(k+3)-vital_sign(k-3)))/(32*Tf);
end

for k=4:frame-6
   vital_sign_diff2(k-3)=(4*vital_sign(k)+(vital_sign(k+1)+vital_sign(k-1))-2*(vital_sign(k+2)+vital_sign(k-2))-(vital_sign(k+3)+vital_sign(k-3)))/(16*Tf*Tf);
end

%ͨ������Ƶ��0.9-2Hz��ͨ�˲���
filter_heart_diff1=filter(bpf_heart,vital_sign_diff1);
filter_heart_diff2=filter(bpf_heart,vital_sign_diff2);

%��ԭʼ�ź���fft
vital_sign_fft=fft(vital_sign,N);

%˫�ߴ��ź�תΪ���ߴ�
P2 = abs(vital_sign_fft/(N-1));
P1 = P2(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1(2:end-1) = 2*P1(2:end-1);

%�ж�ʹ��һ��΢�ֻ��Ƕ���΢��
threshold=6; %��ֵ
if max(P1)>threshold*mean(P1)
    %�����ź�ʱ��ͼ
    index=1:1:frame-9;
    index=index*Tf;
    figure(1);
    plot(index,filter_heart_diff1);
    xlabel('Time(s)');
    title('�����ź�ʱ��ͼ');

    %�����
    heart_xcorr=xcorr(filter_heart_diff1);
    
    heart_fft=fft(heart_xcorr,N);
    freq=(0:1:N/2)/Tf/N;  %���������źŲ�����Ϊ֡��
    
    %˫�ߴ���Ϊ���ߴ�
    P2_heart = abs(heart_fft/(N-1));
    P1_heart = P2_heart(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
    P1_heart(2:end-1) = 2*P1_heart(2:end-1);
    
    %�����ź�Ƶ��ͼ
    figure(2);
    plot(freq,P1_heart);
    xlim([0,4]);
    xlabel('Frequency(Hz)');
    title('�����ź�Ƶ��ͼ');
else
    %�����ź�ʱ��ͼ
    index=1:1:frame-9;
    index=index*Tf;
    figure(1);
    plot(index,filter_heart_diff2);
    xlabel('Time(s)');
    title('�����ź�ʱ��ͼ');

    %�����
    heart_xcorr=xcorr(filter_heart_diff2);
    
    heart_fft=fft(heart_xcorr,N);
    freq=(0:1:N/2)/Tf/N;  %���������źŲ�����Ϊ֡��
    
    %˫�ߴ���Ϊ���ߴ�
    P2_heart = abs(heart_fft/(N-1));
    P1_heart = P2_heart(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
    P1_heart(2:end-1) = 2*P1_heart(2:end-1);
    
    %�����ź�Ƶ��ͼ
    figure(2);
    plot(freq,P1_heart);
    xlim([0,4]);
    xlabel('Frequency(Hz)');
    title('�����ź�Ƶ��ͼ');
end
