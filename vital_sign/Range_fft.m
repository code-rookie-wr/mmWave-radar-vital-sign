clear;
[retVal] = readDCA1000_1('./demo.bin');
global numChirps;
global numADCSamples;%��������
RX1data = reshape(retVal(1,:),numADCSamples,numChirps);   %RX1����
RX2data = reshape(retVal(2,:),numADCSamples,numChirps);   %RX2
RX3data = reshape(retVal(3,:),numADCSamples,numChirps);   %RX3
RX4data = reshape(retVal(4,:),numADCSamples,numChirps);   %RX4

c=3.0e8;  
slope=60e12; %��Ƶб��
Tc=50e-6;      %chirp����
B=slope*Tc;    %��Ƶ����
Fs=4e6;        %������
f0=60.36e9;       %��ʼƵ��
lambda=c/f0;   %�״��źŲ���
d=lambda/2;    %�������м��
frame=400;     %֡��
Tf=0.05;       %֡����

data = RX1data(:,(1:2:end)); %ȡ��һ���еĵ�1��256�е�����(һ��chirp)
index = 1:1:numADCSamples; 

%���ɴ�
range_win = hamming(numADCSamples);
%�Ӵ�����
din_win = data.* range_win;
%fft����
datafft = fft(din_win);
%samplesת��Ϊfreq
freq_bin = (index-1) * Fs / numADCSamples;
%freqת��Ϊrange
range_bin = freq_bin * c / 2 / slope;
dBFs=20*log10(abs(datafft)/2^16);
%ʱ��ָ��
Time=Tf*frame:-Tf:Tf;
%��������ͼ
figure(1);
imagesc(range_bin,Time,abs(datafft)');
colorbar;
%xlim([0,5]);
% title('ʱ�����ͼ');
xlabel('Distance(meters)');
ylabel('time(s)');

