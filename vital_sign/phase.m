clear;
[retVal] = readDCA1000_1('C:\Users\Administrator\Desktop\�״�����\0.3-1wr.bin');
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
f0=60.36e9;       %��ʼƵ��
lambda=c/f0;   %�״��źŲ���
d=lambda/2;    %�������м��
frame=400;     %֡��
Tf=0.05;       %֡����

range_win = hamming(numADCSamples); %���ɺ�����
for k=1:1:frame
    din_win(:,k)=RX1data(:,2*k-1).*range_win; %���ź���Range-fft
    datafft(:,k)=fft(din_win(:,k));
end
%�ҵ�Range-bin��ֵ
rangeBinStartIndex=8;
rangeBinEndIndex=24;
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
thresh=1.5;
for k=1:frame-3
    phaseUsedComputation(:,k)=filter_RemoveImpulseNoise(delta_phase(:,k),delta_phase(:,k+1),delta_phase(:,k+2),thresh);
end
index=1:1:frame-3;
index=index*Tf;
%ԭʼ�źŴ�ͨ�˲�
filter_delta_phase=filter(bpf_vitalsign,phaseUsedComputation);
%����λת��Ϊλ��
% delta_d=(filter_delta_phase*lambda)/4/pi;
figure(1);
%ԭʼ�ź���λʱ��ͼ
plot(index,filter_delta_phase);
xlabel('Time(s)');
ylabel('Amplitude');
title('ԭʼ�ź�ʱ��ͼ');
%��ԭʼ��λ�ź���fft
filter_delta_phase_fft=fft(filter_delta_phase);
%˫�ߴ��ź�תΪ���ߴ�
freq=(0:1:frame/2)/Tf/frame;  %���������źŲ�����Ϊ֡��
P2 = abs(filter_delta_phase_fft/(frame-1));
P1 = P2(1:frame/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1(2:end-1) = 2*P1(2:end-1);
%ԭʼ�ź���λƵ��ͼ
figure(2);
plot(freq,P1);
xlim([0,2]);
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('ԭʼ�ź�Ƶ��ͼ');
%�����źŴ�ͨ�˲�
filter_delta_phase_breathe=filter(bpf_breathe,phaseUsedComputation);
%�����ź���λʱ��ͼ
figure(3);
plot(index,filter_delta_phase_breathe);
xlabel('Time(s)');
ylabel('Amplitude');
title('�����ź�ʱ��ͼ');
%��ԭʼ��λ�ź���fft
filter_delta_phase_breathe_fft=fft(filter_delta_phase_breathe);
%˫�ߴ��ź�תΪ���ߴ�
P2_breathe = abs(filter_delta_phase_breathe_fft/(frame-1));
P1_breathe = P2_breathe(1:frame/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1_breathe(2:end-1) = 2*P1_breathe(2:end-1);
%�����ź���λƵ��ͼ
figure(4);
plot(freq,P1_breathe);
xlim([0,2]);
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('�����ź�Ƶ��ͼ');
%�����źŴ�ͨ�˲�
filter_delta_phase_heart=filter(bpf_heart,phaseUsedComputation);
%�����ź���λʱ��ͼ
figure(5);
plot(index,filter_delta_phase_heart);
xlabel('Time(s)');
ylabel('Amplitude');
title('�����ź�ʱ��ͼ');
%��ԭʼ��λ�ź���fft
delta_phase_heart_fft=fft(filter_delta_phase_heart);
%˫�ߴ��ź�תΪ���ߴ�
P2_heart = abs(delta_phase_heart_fft/(frame-1));
P1_heart = P2_heart(1:frame/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1_heart(2:end-1) = 2*P1_heart(2:end-1);
%�����ź���λƵ��ͼ
figure(6);
plot(freq,P1_heart);
xlim([0,4]);
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('�����ź�Ƶ��ͼ');
%emd
% figure(7);
% imf=eemd(filter_delta_phase,0.4,50);
% [a,b]=size(imf);
% for i=1:b
% %     suptitle('Empirical Mode Decomposition');
%     subplot(b,2,2*i-1);
%     plot(index,imf(:,i));
%     xlabel('Time(s)');
%     subplot(b,2,2*i);
%     imf_fft=fft(imf(:,i));
%     P2_imf = abs(imf_fft/(frame-1));
%     P1_imf = P2_imf(1:frame/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
%     P1_imf(2:end-1) = 2*P1_imf(2:end-1);
%     plot(freq,P1_imf);
%     xlim([0,2]);
%     xlabel('Frequency(Hz)');
% end