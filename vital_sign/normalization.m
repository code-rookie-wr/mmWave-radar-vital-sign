function [Output] = normalization(Input, lowMargin, highMargin)
%Input:����
%lowMargin:��һ���±߽�
%highMargin:��һ���ϱ߽�
%Output:���
maxValue = max(Input);
minValue = min(Input);
Output = (Input-minValue)/(maxValue-minValue);
Output = Output*(highMargin-lowMargin)+lowMargin;

end