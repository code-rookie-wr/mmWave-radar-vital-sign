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

%真实胸腔位移信号
% for i=1:frame-3
%    vital_sign(i)= (phaseUsedComputation(1)+phaseUsedComputation(i))*i/2;
% end

vital_sign=phaseUsedComputation;
vital_sign=filter(bpf_vitalsign,vital_sign);

%VMD参数设置
alpha = 5000;      % moderate bandwidth constraint
tau = 0;            % noise-tolerance (no strict fidelity enforcement)0
K = 3;              % modes
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-6;

%VMD
[u, u_hat, omega] = VMD(vital_sign, alpha, tau, K, DC, init, tol);

%% 呼吸信号
%呼吸信号时域图
respiration=filter(bpf_breathe,u(1,:));
index=1:1:frame-3;
index=index*Tf;
figure(3);
plot(index,respiration,'b');
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('呼吸信号','FontWeight','bold');

%呼吸信号频域图
freq=(0:1:N/2)/Tf/N;
figure(4);
imf_breathe_fft=fft(respiration,N);
P2_breathe_imf = abs(imf_breathe_fft/(N-1));
P1_breathe_imf = P2_breathe_imf(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1_breathe_imf(2:end-1) = 2*P1_breathe_imf(2:end-1);
plot(freq,P1_breathe_imf);
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('呼吸信号频谱图','FontWeight','bold');

%% 呼吸气面积比
[inhale_to_exhale, exhale_to_inhale, point30Percent, point70Percent] = exhale_inhale_area_ratio(respiration,Tf);

%将30%-70%的点连线
figure(3);
hold on;
plot(point30Percent(1, :), point30Percent(2, :), 'r-x');
plot(point70Percent(1, :), point70Percent(2, :), 'g-o');
hold off;

%% 心跳信号
%心跳信号时域图
heart=filter(bpf_heart,u(2,:));
figure(5);
plot(index,heart);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold ');
title('心跳信号','FontWeight','bold');

%心跳信号频域图
figure(6);
imf_heart_fft=fft(heart,N);
P2_heart_imf = abs(imf_heart_fft/(N-1));
P1_heart_imf = P2_heart_imf(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1_heart_imf(2:end-1) = 2*P1_heart_imf(2:end-1);
plot(freq,P1_heart_imf);
xlim([0,4]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('心跳信号频谱图','FontWeight','bold');

%% HRV
% [IBI, MEAN, SDNN, r_MSSD] = mmHRV(heart, Tf);
% figure(7);
% plot(IBI(1,:));
% axis auto;
% xlabel('IBI Index','FontWeight','bold');
% ylabel('IBI(s)','FontWeight','bold');
% title('mmHRV的IBI估计','FontWeight','bold');

%% 原始信号
%生命体征信号重构
vital_sign=filter(bpf_heart,u(2,:))+filter(bpf_breathe,u(1,:));

%生命体征信号时域图
figure(1);
plot(index,vital_sign);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('心肺信号','FontWeight','bold');

%生命体征信号频域图
figure(2);
imf_vital_sign_fft=fft(vital_sign,N);
P2_vital_sign_imf = abs(imf_vital_sign_fft/(N-1));
P1_vital_sign_imf = P2_vital_sign_imf(1:N/2+1);   %此时选取前半部分，因为fft之后为对称的双边谱
P1_vital_sign_imf(2:end-1) = 2*P1_vital_sign_imf(2:end-1);
plot(freq,P1_vital_sign_imf);
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('心肺信号频谱图','FontWeight','bold');
