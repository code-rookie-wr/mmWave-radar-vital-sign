function [Output] = normalization(Input, lowMargin, highMargin)
%Input:输入
%lowMargin:归一化下边界
%highMargin:归一化上边界
%Output:输出
maxValue = max(Input);
minValue = min(Input);
Output = (Input-minValue)/(maxValue-minValue);
Output = Output*(highMargin-lowMargin)+lowMargin;

end