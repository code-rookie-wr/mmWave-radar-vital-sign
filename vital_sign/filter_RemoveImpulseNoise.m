function y = filter_RemoveImpulseNoise(dataPrev2,dataPrev1,dataCurr,thresh)
pDataIn=[];
pDataIn(1)=dataPrev2;
pDataIn(2)=dataPrev1;
pDataIn(3)=dataCurr;

backwardDiff=pDataIn(2)-pDataIn(1);
forwardDiff=pDataIn(2)-pDataIn(3);

x1=0;
x2=2;
y1=pDataIn(1);
y2=pDataIn(3);
x=1;

if( ((forwardDiff > thresh) && (backwardDiff > thresh)) || ((forwardDiff < -thresh) && (backwardDiff < -thresh)) )
    y=y1+( ((x-x1)*(y2-y1))/(x2-x1) );
else
    y=pDataIn(2);
end
end

