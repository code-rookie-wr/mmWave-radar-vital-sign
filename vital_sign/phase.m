% 该文件仅用于对信号相位变化的监测，用于为进一步算法做铺垫
clear;
[retVal] = readDCA1000_1('./demo.bin');
global numChirps;
global numADCSamples;%采样点数
RX1data = reshape(retVal(1,:),numADCSamples,numChirps);   %RX1数据
RX2data = reshape(retVal(2,:),numADCSamples,numChirps);   %RX2
RX3data = reshape(retVal(3,:),numADCSamples,numChirps);   %RX3
RX4data = reshape(retVal(4,:),numADCSamples,numChirps);   %RX4

c=3.0e8;  
slope=60e12; %调频斜率
Tc=50e-6;      %chirp周期
B=slope*Tc;    %调频带宽
Fs=4e6;        %采样率
f0=60.36e9;       %初始频率
lambda=c/f0;   %雷达信号波长
d=lambda/2;    %天线阵列间距
frame=400;     %帧数
Tf=0.05;       %帧周期

range_win = hamming(numADCSamples); %生成海明窗
for k=1:1:frame
    din_win(:,k)=RX1data(:,2*k-1).*range_win; %对信号做Range-fft
    datafft(:,k)=fft(din_win(:,k));
end
%找到Range-bin峰值
rangeBinStartIndex=8;
rangeBinEndIndex=24;
for k=1:1:frame
    for j=rangeBinStartIndex:1:numADCSamples 
        if(abs(datafft(j,k))==max(abs(datafft((rangeBinStartIndex:rangeBinEndIndex),k)))) 
            data(:,k)=datafft(j,k);
        end
    end
end
%获取信号实部、虚部
for k=1:frame
    data_real(:,k)=real(data(:,k));
    data_imag(:,k)=imag(data(:,k));
end
%计算信号相位
for k=1:frame
    signal_phase(:,k)=atan(data_imag(:,k)/data_real(:,k));
end
%相位展开
for k=2:frame
    diff=signal_phase(:,k)-signal_phase(:,k-1);
    if diff>pi/2
        signal_phase(:,(k:end))=signal_phase(:,(k:end))-pi;
    elseif diff<-pi/2
        signal_phase(:,(k:end))=signal_phase(:,(k:end))+pi;
    end
end
%计算相位差
for k=1:frame-1
    delta_phase(:,k)=signal_phase(:,k+1)-signal_phase(:,k);
end
%从波形中去除脉冲噪声
thresh=1.5;
for k=1:frame-3
    phaseUsedComputation(:,k)=filter_RemoveImpulseNoise(delta_phase(:,k),delta_phase(:,k+1),delta_phase(:,k+2),thresh);
end
index=1:1:frame-3;
index=index*Tf;
%原始信号带通滤波
filter_delta_phase=filter(bpf_vitalsign,phaseUsedComputation);
%将相位转换为位移
% delta_d=(filter_delta_phase*lambda)/4/pi;
figure(1);
%原始信号相位时域图
plot(index,filter_delta_phase);
xlabel('Time(s)');
ylabel('Amplitude');
title('原始信号时域图');
%对原始相位信号做fft
filter_delta_phase_fft=fft(filter_delta_phase);
%双边带信号转为单边带
freq=(0:1:frame/2)/Tf/frame;  %生命特征信号采样率为帧数
P2 = abs(filter_delta_phase_fft/(frame-1));
P1 = P2(1:frame/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1(2:end-1) = 2*P1(2:end-1);
%原始信号相位频域图
figure(2);
plot(freq,P1);
xlim([0,2]);
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('原始信号频域图');
%呼吸信号带通滤波
filter_delta_phase_breathe=filter(bpf_breathe,phaseUsedComputation);
%呼吸信号相位时域图
figure(3);
plot(index,filter_delta_phase_breathe);
xlabel('Time(s)');
ylabel('Amplitude');
title('呼吸信号时域图');
%对原始相位信号做fft
filter_delta_phase_breathe_fft=fft(filter_delta_phase_breathe);
%双边带信号转为单边带
P2_breathe = abs(filter_delta_phase_breathe_fft/(frame-1));
P1_breathe = P2_breathe(1:frame/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1_breathe(2:end-1) = 2*P1_breathe(2:end-1);
%呼吸信号相位频域图
figure(4);
plot(freq,P1_breathe);
xlim([0,2]);
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('呼吸信号频域图');
%心跳信号带通滤波
filter_delta_phase_heart=filter(bpf_heart,phaseUsedComputation);
%心跳信号相位时域图
figure(5);
plot(index,filter_delta_phase_heart);
xlabel('Time(s)');
ylabel('Amplitude');
title('心跳信号时域图');
%对原始相位信号做fft
delta_phase_heart_fft=fft(filter_delta_phase_heart);
%双边带信号转为单边带
P2_heart = abs(delta_phase_heart_fft/(frame-1));
P1_heart = P2_heart(1:frame/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1_heart(2:end-1) = 2*P1_heart(2:end-1);
%心跳信号相位频域图
figure(6);
plot(freq,P1_heart);
xlim([0,4]);
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('心跳信号频域图');
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
%     P1_imf = P2_imf(1:frame/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
%     P1_imf(2:end-1) = 2*P1_imf(2:end-1);
%     plot(freq,P1_imf);
%     xlim([0,2]);
%     xlabel('Frequency(Hz)');
% end