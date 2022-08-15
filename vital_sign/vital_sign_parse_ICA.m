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
index=1:1:frame-3;
index=index*Tf;

%fast-ICA
[z]=fastICA(phaseUsedComputation,1);
nmf_result=z;

%原始信号带通滤波
filter_delta_phase=filter(bpf_vitalsign,nmf_result);

%将相位转换为位移
% delta_d=(filter_delta_phase*lambda)/4/pi;

%真实胸腔位移信号
% for i=1:frame-3
%    vital_sign(i)= (filter_delta_phase(1)+filter_delta_phase(i))*i/2;
% end

vital_sign=filter_delta_phase;

%原始信号时域图
figure(1);
plot(index,vital_sign);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('心肺信号','FontWeight','bold');

%对原始信号做fft
vital_sign_fft=fft(vital_sign,N);

%双边带信号转为单边带
freq=(0:1:N/2)/Tf/N;  %生命特征信号采样率为帧数
P2 = abs(vital_sign_fft/(N-1));
P1 = P2(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1(2:end-1) = 2*P1(2:end-1);

%原始信号频域图
figure(2);
plot(freq,P1);
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('心肺信号频谱图','FontWeight','bold');

%呼吸信号带通滤波
filter_delta_phase_breathe=filter(bpf_breathe,vital_sign);

%信号转为真实呼吸信号
% for i=1:frame-3
%    breathe(i)= (filter_delta_phase_breathe(1)+filter_delta_phase_breathe(i))*i/2;
% end

breathe=filter_delta_phase_breathe;

%呼吸信号时域图
figure(3);
plot(index,breathe);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('呼吸信号','FontWeight','bold');

%对呼吸信号做fft
breathe_fft=fft(breathe,N);

%双边带信号转为单边带
P2_breathe = abs(breathe_fft/(N-1));
P1_breathe = P2_breathe(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1_breathe(2:end-1) = 2*P1_breathe(2:end-1);

%呼吸信号频域图
figure(4);
plot(freq,P1_breathe);
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('呼吸信号频谱图','FontWeight','bold');


%心跳信号带通滤波
filter_delta_phase_heart=filter(bpf_heart,vital_sign);

%真实心跳位移信号
% for i=1:frame-3
%    heart(i)= (filter_delta_phase_heart(1)+filter_delta_phase_heart(i))*i/2;
% end

heart=filter_delta_phase_heart;

%心跳信号时域图
figure(5);
plot(index,heart);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold ');
title('心跳信号','FontWeight','bold');

%对心跳信号做fft
heart_fft=fft(heart,N);

%双边带信号转为单边带
P2_heart = abs(heart_fft/(N-1));
P1_heart = P2_heart(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1_heart(2:end-1) = 2*P1_heart(2:end-1);

%心跳谐波检测
[heart_peaks,heart_peaksnum]=findpeaks(P1_heart,0.9,2,N,Tf);%0.9Hz-2Hz
[heart_harmonic_peaks,heart_harmonic_peaksnum]=findpeaks(P1_heart,1.8,4,N,Tf);%1.8Hz-4Hz
heart_peaks=heart_peaks/N/Tf;
heart_harmonic_peaks=heart_harmonic_peaks/N/Tf;
[heart_peaks_row,heart_peaks_column]=size(heart_peaks);
[heart_harmonic_peaks_row,heart_harmonic_peaks_column]=size(heart_harmonic_peaks);
for i=1:heart_peaks_column
    if max(P1_heart)-P1_heart(round(heart_peaks(i)*N*Tf)+1)<0.3 %当出现两峰值之差在一定范围内，再将找到谐波的信号进行扩大
        for j=1:heart_harmonic_peaks_column
            if heart_harmonic_peaks(j)/heart_peaks(i)==2
                P1_heart(round(heart_peaks(i)*N*Tf)+1)=2*P1_heart(round(heart_peaks(i)*N*Tf)+1);
            end
        end
    end
end

%心跳信号频域图
figure(6);
plot(freq,P1_heart);
xlim([0,4]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('心跳信号频谱图','FontWeight','bold');
