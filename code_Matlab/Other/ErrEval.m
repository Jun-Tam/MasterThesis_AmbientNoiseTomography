function [RdErrStd,RdErrMean]=...
    ErrEval(VphPthMdlIn,StPos,StDist,StAzim,DistRng,AzimRng,NumThr)

tic;
NumFreq=size(VphPthMdlIn,2);
RdErrStd=zeros(500,NumFreq);
RdErrMean=zeros(1,NumFreq);
for i=1:NumFreq
    NumPair=find(VphPthMdlIn(:,i)~=0, 1, 'last' );
    VelStr=zeros(size(VphPthMdlIn,1),50);
    for j=1:NumPair
        n=1;
        LonGr1=mean(StPos(j,[1 3],i));
        LatGr1=mean(StPos(j,[2 4],i));
        for k=1:NumPair
            LonGr2=mean(StPos(k,[1 3],i));
            LatGr2=mean(StPos(k,[2 4],i));            
            [arclenGr,~]=distance(LatGr1,LonGr1,LatGr2,LonGr2);            
            DistGr=sqrt((arclenGr*2*pi/360*earthRadius)^2)/1000;
            DistDif=abs(StDist(j,i)-StDist(k,i));
            AzimDif=abs(StAzim(j,i)-StAzim(k,i));      
            if DistGr + DistDif < DistRng && AzimDif < AzimRng
                VelStr(j,n)=VphPthMdlIn(k,i);
                n=n+1;
            end
        end
    end
    m=1;
    for j=1:NumPair
        IdxNon0=VelStr(j,:)~=0;
        if sum(IdxNon0) > NumThr
            RdErrStd(m,i)=std(VelStr(j,IdxNon0));
            m=m+1;
        end
    end
    RdErrMean(i)=mean(RdErrStd(1:m-1,i));
end
LapseTime=toc;
display(strcat('Error Evaluation:',num2str(LapseTime),'sec'));

end
