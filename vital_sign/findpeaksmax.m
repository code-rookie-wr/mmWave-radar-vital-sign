function peaks=findpeaksmax(data,StartFreq,EndFreq,fft_size,Tf)

StartIndex=round(StartFreq*fft_size*Tf);
EndIndex=round(EndFreq*fft_size*Tf);
for i=StartIndex:EndIndex
    if (data(i)==max(data))
       peaks=i-1;
    end
end