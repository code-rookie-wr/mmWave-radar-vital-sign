function [index, xt] = go_cfar(xc, N, pro_N, PAD)

% N表示参考窗口的长度
% pro_N表示用于检测门限的百分比
% PAD用于抵消概率密度函数中的偏差

alpha=N.*(PAD.^(-1./N)-1);
index=1+N/2+pro_N/2:length(xc)-N/2-pro_N/2;
XT=zeros(1,length(index));

for i=index
    cell_left=xc(1,i-N/2-pro_N/2:i-pro_N/2-1);
    cell_right=xc(1,i+pro_N/2+1:i+N/2+pro_N/2);
    Z=max([mean(cell_left),mean(cell_right)]);

    xt(1,i-N/2-pro_N/2)=Z.*alpha;
end

end