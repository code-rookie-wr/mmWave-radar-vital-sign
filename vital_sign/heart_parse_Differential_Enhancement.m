clear;
[retVal] = readDCA1000_1('./demo.bin');
global numChirps;
global numADCSamples;%采样点数
RX1data = reshape(retVal(1,:),numADCSamples,numChirps);   %RX1数据
RX2data = reshape(retVal(2,:),numADCSamples,numChirps);   %RX2
RX3data = reshape(retVal(3,:),numADCSamples,numChirps);   %RX3
RX4data = reshape(retVal(4,:),numADCSamples,numChirps);   %RX4

c=3.0e8;  
slope=59.997e12; %调频斜率
Tc=50e-6;      %chirp周期
B=slope*Tc;    %调频带宽
Fs=4e6;        %采样率
f0=60.36e9;    %初始频率
lambda=c/f0;   %雷达信号波长
d=lambda/2;    %天线阵列间距
frame=400;     %帧数
Tf=0.05;       %帧周期
N=1024;        %FFT点数

%加海明窗
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
thresh=0.8;
for k=1:frame-3
    phaseUsedComputation(:,k)=filter_RemoveImpulseNoise(delta_phase(:,k),delta_phase(:,k+1),delta_phase(:,k+2),thresh);
end

%原始信号带通滤波
filter_delta_phase=filter(bpf_vitalsign,phaseUsedComputation);

vital_sign=filter_delta_phase;

%计算差分时间序列
for k=4:frame-6
   vital_sign_diff1(k-3)=(5*(vital_sign(k+1)-vital_sign(k-1))+4*(vital_sign(k+2)-vital_sign(k-2))+(vital_sign(k+3)-vital_sign(k-3)))/(32*Tf);
end

for k=4:frame-6
   vital_sign_diff2(k-3)=(4*vital_sign(k)+(vital_sign(k+1)+vital_sign(k-1))-2*(vital_sign(k+2)+vital_sign(k-2))-(vital_sign(k+3)+vital_sign(k-3)))/(16*Tf*Tf);
end

%通过心跳频段0.9-2Hz带通滤波器
filter_heart_diff1=filter(bpf_heart,vital_sign_diff1);
filter_heart_diff2=filter(bpf_heart,vital_sign_diff2);

%对原始信号做fft
vital_sign_fft=fft(vital_sign,N);

%双边带信号转为单边带
P2 = abs(vital_sign_fft/(N-1));
P1 = P2(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1(2:end-1) = 2*P1(2:end-1);

%判定使用一阶微分还是二阶微分
threshold=6; %阈值
if max(P1)>threshold*mean(P1)
    %心跳信号时域图
    index=1:1:frame-9;
    index=index*Tf;
    figure(1);
    plot(index,filter_heart_diff1);
    xlabel('Time(s)');
    title('心跳信号时域图');

    %自相关
    heart_xcorr=xcorr(filter_heart_diff1);
    
    heart_fft=fft(heart_xcorr,N);
    freq=(0:1:N/2)/Tf/N;  %生命特征信号采样率为帧数
    
    %双边带变为单边带
    P2_heart = abs(heart_fft/(N-1));
    P1_heart = P2_heart(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
    P1_heart(2:end-1) = 2*P1_heart(2:end-1);
    
    %心跳信号频域图
    figure(2);
    plot(freq,P1_heart);
    xlim([0,4]);
    xlabel('Frequency(Hz)');
    title('心跳信号频域图');
else
    %心跳信号时域图
    index=1:1:frame-9;
    index=index*Tf;
    figure(1);
    plot(index,filter_heart_diff2);
    xlabel('Time(s)');
    title('心跳信号时域图');

    %自相关
    heart_xcorr=xcorr(filter_heart_diff2);
    
    heart_fft=fft(heart_xcorr,N);
    freq=(0:1:N/2)/Tf/N;  %生命特征信号采样率为帧数
    
    %双边带变为单边带
    P2_heart = abs(heart_fft/(N-1));
    P1_heart = P2_heart(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
    P1_heart(2:end-1) = 2*P1_heart(2:end-1);
    
    %心跳信号频域图
    figure(2);
    plot(freq,P1_heart);
    xlim([0,4]);
    xlabel('Frequency(Hz)');
    title('心跳信号频域图');
end
