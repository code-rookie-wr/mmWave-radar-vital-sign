function [IBI, MEAN, SDNN, r_MSSD] = mmHRV(heart, Tf)

%获取极值点
extrMaxIndex = find(diff(sign(diff(heart)))==-2)+1;
extrMinIndex = find(diff(sign(diff(heart)))==+2)+1;

%计算IBI,第一行表示极大值算出的结果，第二行表示极小值算出的结果
[m, n] = size(extrMaxIndex);
for i = 1:n-1
    IBI(1,i) = (extrMaxIndex(i+1)-extrMaxIndex(i))*Tf;
end
[m, n] = size(extrMinIndex);
for i = 1:n-1
    IBI(2,i) = (extrMinIndex(i+1)-extrMinIndex(i))*Tf;
end

%计算MEAN
MEAN = mean(IBI, 2);

%计算SDNN
SDNN = std(IBI, 1, 2);

%计算r_MSSD
[m, n] = size(IBI);
for i = 1:n-1
    delta_RR(:,i) = IBI(:,i+1)-IBI(:,i);
end
r_MSSD = sqrt(sum(delta_RR,2).^2/(n-1));

end