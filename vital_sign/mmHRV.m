function [IBI, MEAN, SDNN, r_MSSD] = mmHRV(heart, Tf)

%��ȡ��ֵ��
extrMaxIndex = find(diff(sign(diff(heart)))==-2)+1;
extrMinIndex = find(diff(sign(diff(heart)))==+2)+1;

%����IBI,��һ�б�ʾ����ֵ����Ľ�����ڶ��б�ʾ��Сֵ����Ľ��
[m, n] = size(extrMaxIndex);
for i = 1:n-1
    IBI(1,i) = (extrMaxIndex(i+1)-extrMaxIndex(i))*Tf;
end
[m, n] = size(extrMinIndex);
for i = 1:n-1
    IBI(2,i) = (extrMinIndex(i+1)-extrMinIndex(i))*Tf;
end

%����MEAN
MEAN = mean(IBI, 2);

%����SDNN
SDNN = std(IBI, 1, 2);

%����r_MSSD
[m, n] = size(IBI);
for i = 1:n-1
    delta_RR(:,i) = IBI(:,i+1)-IBI(:,i);
end
r_MSSD = sqrt(sum(delta_RR,2).^2/(n-1));

end