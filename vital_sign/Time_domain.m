% ���ļ�������I/Q��·���ź��Լ�ȡģ����ź�
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
f0=60.36e9;       %��ʼƵ��
lambda=c/f0;   %�״��źŲ���
d=lambda/2;    %�������м��
frame=500;     %֡��
Tf=0.05;       %֡����

%ʱ��ͼ
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
