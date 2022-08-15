% 此文件用于对于长时间呼吸速率以及心率的检测，此demo仅有20s数据，这里仅供参考
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
f0=60.36e9;    %初始频率
lambda=c/f0;   %雷达信号波长
d=lambda/2;    %天线阵列间距
frame=400;     %帧数
Tf=0.05;       %帧周期
N=1024;        %FFT点数

range_win = hamming(numADCSamples); %生成海明窗
for k=1:1:frame
    din_win(:,k)=RX1data(:,2*k-1).*range_win; %对信号做Range-fft
    datafft(:,k)=fft(din_win(:,k));
end

%找到Range-bin峰值
rangeBinStartIndex=3;%分辨率0.1m
rangeBinEndIndex=10;
for k=1:1:frame
    for j=rangeBinStartIndex:1:numADCSamples 
        if(abs(datafft(j,k))==max(abs(datafft((rangeBinStartIndex:rangeBinEndIndex),k)))) 
            data(:,k)=datafft(j,k);
        end
    end
end

%对数据加个300帧窗口，每次以20帧进行滑动
frame_process=300;
slide_width=20;
for i=0:((frame-frame_process)/slide_width)
    data_process(i+1,:)=data(:,(1+slide_width*i):(frame_process+slide_width*i));

    
%获取信号实部、虚部
for k=1:frame_process
    data_real(i+1,k)=real(data_process(i+1,k));
    data_imag(i+1,k)=imag(data_process(i+1,k));
end

%计算信号相位
for k=1:frame_process
    signal_phase(i+1,k)=atan(data_imag(i+1,k)/data_real(i+1,k));
end

%相位展开
for k=2:frame_process
    diff=signal_phase(i+1,k)-signal_phase(i+1,k-1);
    if diff>pi/2
        signal_phase(i+1,(k:end))=signal_phase(i+1,(k:end))-pi;
    elseif diff<-pi/2
        signal_phase(i+1,(k:end))=signal_phase(i+1,(k:end))+pi;
    end
end

%计算相位差
for k=1:frame_process-1
    delta_phase(i+1,k)=signal_phase(i+1,k+1)-signal_phase(i+1,k);
end

%从波形中去除脉冲噪声
thresh=1;
for k=1:frame_process-3
    phaseUsedComputation(i+1,k)=filter_RemoveImpulseNoise(delta_phase(i+1,k),delta_phase(i+1,k+1),delta_phase(i+1,k+2),thresh);
end

%呼吸信号带通滤波
filter_delta_phase_breathe(i+1,:)=filter(bpf_breathe,phaseUsedComputation(i+1,:));

%信号转为真实呼吸信号
% for k=1:frame_process-3
%    breathe(i+1,k)= (filter_delta_phase_breathe(1)+filter_delta_phase_breathe(i+1,k))*k/2;
% end

breathe=filter_delta_phase_breathe;

%对呼吸信号做fft
breathe_fft(i+1,:)=fft(breathe(i+1,:),N);

%双边带信号转为单边带
P2_breathe = abs(breathe_fft(i+1,:)/(N-1));
P1_breathe = P2_breathe(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1_breathe(2:end-1) = 2*P1_breathe(2:end-1);
ssb_breathe(i+1,:)=P1_breathe;

%找到呼吸频域峰值
breathe_max(i+1)=findpeaksmax(ssb_breathe(i+1,:),0.1,0.6,N,Tf)/N/Tf*60;

%心跳信号带通滤波
filter_delta_phase_heart(i+1,:)=filter(bpf_heart,phaseUsedComputation(i+1,:));

%真实心跳位移信号
% for k=1:frame_process-3
%    heart(i+1,k)= (filter_delta_phase_heart(1)+filter_delta_phase_heart(i+1,k))*k/2;
% end

heart=filter_delta_phase_heart;

%对心跳信号做fft
heart_fft((i+1),:)=fft(heart((i+1),:),N);

%双边带信号转为单边带
P2_heart = abs(heart_fft((i+1),:)/(N-1));
P1_heart = P2_heart(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1_heart(2:end-1) = 2*P1_heart(2:end-1);

%心跳谐波检测
[heart_peaks,heart_peaksnum]=findpeaks(P1_heart,0.9,2,N,Tf);%0.9Hz-2Hz
[heart_harmonic_peaks,heart_harmonic_peaksnum]=findpeaks(P1_heart,1.8,4,N,Tf);%1.8Hz-4Hz
heart_peaks=heart_peaks/N/Tf;
heart_harmonic_peaks=heart_harmonic_peaks/N/Tf;
[heart_peaks_row,heart_peaks_column]=size(heart_peaks);
[heart_harmonic_peaks_row,heart_harmonic_peaks_column]=size(heart_harmonic_peaks);
for k=1:heart_peaks_column
    if max(P1_heart)-P1_heart(round(heart_peaks(k)*N*Tf)+1)<0.3 %当出现两峰值之差在一定范围内，再将找到谐波的信号进行扩大
        for j=1:heart_harmonic_peaks_column
            if heart_harmonic_peaks(j)/heart_peaks(k)==2
                P1_heart(round(heart_peaks(k)*N*Tf)+1)=2*P1_heart(round(heart_peaks(k)*N*Tf)+1);
            end
        end
    end
end
ssb_heart(i+1,:)=P1_heart;

%找到心跳频域峰值
heart_max(i+1)=findpeaksmax(ssb_heart(i+1,:),0.9,2,N,Tf)/N/Tf*60;
end

%设定时间指数
index=(frame_process/20):(slide_width/20):(frame/20);

%绘制呼吸动态变化图
figure(1);
plot(index,breathe_max);
xlabel('time(s)');
ylabel('times/min');
title('呼吸次数动态变化图');

%绘制心跳动态变化图
figure(2);
plot(index,heart_max);
xlabel('time(s)');
ylabel('times/min');
title('心跳次数动态变化图');