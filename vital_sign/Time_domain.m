% 该文件绘制了I/Q两路的信号以及取模后的信号
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
f0=60.36e9;       %初始频率
lambda=c/f0;   %雷达信号波长
d=lambda/2;    %天线阵列间距
frame=500;     %帧数
Tf=0.05;       %帧周期

%时域图
figure(1);
index=1:1:numADCSamples;
index=index/1e7;
[m,n]=size(RX1data);
len = m;
plot(index,real(RX1data(:,1)),'b-');hold on;
plot(index,imag(RX1data(:,1)),'r-');
grid on;
legend('real','imag');
figure(2);
plot(index,abs(RX1data(:,1)));
xlim([0 len/1e7]);
%set(gca,'xticklabel',[0:0.5:len]);
xlabel('Time(seconds)');
ylabel('Codes');
title('Time domain plot');
energy=sum(abs(RX1data(:,1)).^2);
