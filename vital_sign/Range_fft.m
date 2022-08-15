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

data = RX1data(:,(1:2:end)); %取第一列中的第1至256行的数据(一个chirp)
index = 1:1:numADCSamples; 

%生成窗
range_win = hamming(numADCSamples);
%加窗操作
din_win = data.* range_win;
%fft操作
datafft = fft(din_win);
%samples转换为freq
freq_bin = (index-1) * Fs / numADCSamples;
%freq转换为range
range_bin = freq_bin * c / 2 / slope;
dBFs=20*log10(abs(datafft)/2^16);
%时间指数
Time=Tf*frame:-Tf:Tf;
%绘制热力图
figure(1);
imagesc(range_bin,Time,abs(datafft)');
colorbar;
%xlim([0,5]);
% title('时间距离图');
xlabel('Distance(meters)');
ylabel('time(s)');

