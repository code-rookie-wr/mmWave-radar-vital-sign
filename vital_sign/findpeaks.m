function [peaks,peaksnum]=findpeaks(data,StartFreq,EndFreq,fft_size,Tf)

peaksnum=0;
StartIndex=round(StartFreq*fft_size*Tf);
EndIndex=round(EndFreq*fft_size*Tf);
t=1;
for i=StartIndex:EndIndex
    if (data(i)>data(i-1))&&(data(i)>data(i-2))&&(data(i)>data(i+1))&&(data(i)>data(i+2))
       peaks(t)=i-1;
       t=t+1;
       peaksnum=peaksnum+1;
    end
end