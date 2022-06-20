clear;
[retVal] = readDCA1000_1('E:\�״�����\wr1.bin');
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
f0=60.36e9;    %��ʼƵ��
lambda=c/f0;   %�״��źŲ���
d=lambda/2;    %�������м��
frame=400;     %֡��
Tf=0.05;       %֡����
N=1024;        %FFT����

%�Ӻ�����
range_win = hamming(numADCSamples); %���ɺ�����
for k=1:1:frame
    din_win(:,k)=RX2data(:,2*k-1).*range_win; %���ź���Range-fft
    datafft(:,k)=fft(din_win(:,k));
end

%�ҵ�Range-bin��ֵ
rangeBinStartIndex=3;%�ֱ���0.1m
rangeBinEndIndex=10;
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
thresh=0.8;
for k=1:frame-3
    phaseUsedComputation(:,k)=filter_RemoveImpulseNoise(delta_phase(:,k),delta_phase(:,k+1),delta_phase(:,k+2),thresh);
end

%��ʵ��ǻλ���ź�
% for i=1:frame-3
%    vital_sign(i)= (phaseUsedComputation(1)+phaseUsedComputation(i))*i/2;
% end

vital_sign=phaseUsedComputation;
vital_sign=filter(bpf_vitalsign,vital_sign);

%VMD��������
alpha = 6000;      % moderate bandwidth constraint
tau = 0;            % noise-tolerance (no strict fidelity enforcement)0
K = 3;              % modes
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-6;

%VMD
[u, u_hat, omega] = VMD(vital_sign, alpha, tau, K, DC, init, tol);

%% �����ź�
%�����ź�ʱ��ͼ
respiration=filter(bpf_breathe,u(1,:));
index=1:1:frame-3;
index=index*Tf;
figure(3);
plot(index,respiration,'b');
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('�����ź�','FontWeight','bold');

%�����ź�Ƶ��ͼ
freq=(0:1:N/2)/Tf/N;
figure(4);
imf_breathe_fft=fft(respiration,N);
P2_breathe_imf = abs(imf_breathe_fft/(N-1));
P1_breathe_imf = P2_breathe_imf(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1_breathe_imf(2:end-1) = 2*P1_breathe_imf(2:end-1);
plot(freq,P1_breathe_imf);
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('�����ź�Ƶ��ͼ','FontWeight','bold');

%% �����������
[inhale_to_exhale, exhale_to_inhale, point30Percent, point70Percent] = exhale_inhale_area_ratio(respiration,Tf);

%��30%-70%�ĵ�����
figure(3);
hold on;
plot(point30Percent(1, :), point30Percent(2, :), 'r-x');
plot(point70Percent(1, :), point70Percent(2, :), 'g-o');
hold off;

%% �����ź�
%�����ź�ʱ��ͼ
heart=filter(bpf_heart,u(2,:));
figure(5);
plot(index,heart);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold ');
title('�����ź�','FontWeight','bold');

%�����ź�Ƶ��ͼ
figure(6);
imf_heart_fft=fft(heart,N);
P2_heart_imf = abs(imf_heart_fft/(N-1));
P1_heart_imf = P2_heart_imf(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1_heart_imf(2:end-1) = 2*P1_heart_imf(2:end-1);
plot(freq,P1_heart_imf);
xlim([0,4]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('�����ź�Ƶ��ͼ','FontWeight','bold');

%% HRV
% [IBI, MEAN, SDNN, r_MSSD] = mmHRV(heart, Tf);
% figure(7);
% plot(IBI(1,:));
% axis auto;
% xlabel('IBI Index','FontWeight','bold');
% ylabel('IBI(s)','FontWeight','bold');
% title('mmHRV��IBI����','FontWeight','bold');

%% ԭʼ�ź�
%���������ź��ع�
vital_sign=filter(bpf_heart,u(2,:))+filter(bpf_breathe,u(1,:));

%���������ź�ʱ��ͼ
figure(1);
plot(index,vital_sign,'-k', 'LineWidth', 1.2);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('�ķ��ź�','FontWeight','bold');

%�Ƴ�������߿�
set(gca, 'Visible', 'off');
%���ñ���Ϊ��ɫ
set(gcf, 'color', 'w');

% dataset = (vital_sign(1:end)-min(vital_sign))/(max(vital_sign)-min(vital_sign));
% figure(8);
% for i = 1:318
%     img = plot(index(1:80),dataset(i:(i+79)),'-k');
%     ylim([0,1]);
%     %�Ƴ�������߿�
%     set(gca, 'Visible', 'off');
%     %���ñ���Ϊ��ɫ
%     set(gcf, 'color', 'w');
%     saveas(gcf, ['E:\pythonProject\pyTorch\identity authentication\deep_learning\dataset\wr\', num2str(i+318), '.jpg']);
% end

%���������ź�Ƶ��ͼ
figure(2);
imf_vital_sign_fft=fft(vital_sign,N);
P2_vital_sign_imf = abs(imf_vital_sign_fft/(N-1));
P1_vital_sign_imf = P2_vital_sign_imf(1:N/2+1);   %��ʱѡȡǰ�벿�֣���Ϊfft֮��Ϊ�ԳƵ�˫����
P1_vital_sign_imf(2:end-1) = 2*P1_vital_sign_imf(2:end-1);
plot(freq,P1_vital_sign_imf);
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('�ķ��ź�Ƶ��ͼ','FontWeight','bold');
