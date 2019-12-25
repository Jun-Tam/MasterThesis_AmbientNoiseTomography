function [Flag,VphSyn,VphRef,VrRef,FreqCent,VrVphPth]=...
    DisperCal(Freq0,Vph0,Amp0,VrVph0,Sp0,Dist,FreqRange,FreqStep,...
    FreqBand,PrtbMax,PrtbStep,Sigma,Min0x,Max0x,MinFindPt,PrtbNext,VrThre)

%% Initial setting
% Grid Search
FreqCent=FreqRange(1):FreqStep:FreqRange(2);
PrtbGr=-PrtbMax:PrtbStep:PrtbMax;
NumFreq=length(FreqCent);
NumData=length(Dist);
VphRef=interp1(Freq0,Vph0,FreqCent);
VrRef=interp1(Freq0,VrVph0,FreqCent);
VphSyn=zeros(NumFreq,NumData);
Flag=ones(NumFreq,NumData);
VrVphPth=zeros(NumFreq,NumData);
PrtbFit=zeros(NumFreq,NumData);
% Distance Range
DistBessel=0:0.1:1000;
DistRange=zeros(NumFreq,1);
FreqUse=zeros(NumFreq,NumData);

%% Distance Range (defined by zero crossings of Bessel function)
for i=1:NumFreq
    Bessel=besselj(0,2*pi*FreqCent(i)/VphRef(i)*DistBessel);
    ZeroCross=crossing(Bessel);
    DistRange(i,1)=DistBessel(ZeroCross(Min0x));
    DistRange(i,2)=DistBessel(ZeroCross(Max0x));
end
for i=1:NumData
    if     Dist(i) > max(DistRange(:,1)) && Dist(i) < min(DistRange(:,2))
        FreqUse(:,i)=1;
    elseif Dist(i) > min(DistRange(:,1)) && Dist(i) < max(DistRange(:,1))
        FreqTmp=interp1(DistRange(:,1),FreqCent,Dist(i));
        FreqUse(:,i)=FreqCent > FreqTmp;
    elseif Dist(i) > min(DistRange(:,2)) && Dist(i) < max(DistRange(:,2))
        FreqTmp=interp1(DistRange(:,2),FreqCent,Dist(i));
        FreqUse(:,i)=FreqCent < FreqTmp;
    end
end

%% Dispersion calculation
for i=1:NumData
    jmin=find(FreqUse(:,i)==1, 1 );
    jmax=find(FreqUse(:,i)==1, 1, 'last' );
    if isempty(jmin)~=1
        for j=jmin:jmax
            Idx=Freq0 > FreqCent(j)-FreqBand & Freq0 < FreqCent(j)+FreqBand;
            Data=real((Sp0(i,Idx))');
            Freq=Freq0(Idx);
            Vel=Vph0(Idx);
            Amp=real(Amp0(Idx));
            
            %% Grid Search
            if j==jmin
                Prtb=-PrtbMax:PrtbStep:PrtbMax;
            else
                Prtb=Prtb0-PrtbNext:PrtbStep:Prtb0+PrtbNext;
            end
            [VrVphPth(j,i),PrtbFit(j,i),Prtb0,VphSyn(j,i),FlagTmp]=...
                GrSearch(Prtb,Freq,Vel,Dist(i),Amp,Data,PrtbGr,Idx,Vph0);
            if FlagTmp==0
                break;
            end
        end
    end
end

%% Screening
for i=1:NumData
    for j=1:NumFreq
        if VrVphPth(j,i) < VrThre
            Flag(j,i)=0;
        end
    end
end
for i=1:NumData
    if sum(VphSyn(:,i)~=0) < MinFindPt
        Flag(:,i)=0;
    else
        for j=1:NumFreq
            IdxNonZero=Flag(j,:)==1;
            VphSynMean=mean(VphSyn(j,IdxNonZero));
            Outliar=abs((VphSyn(j,i)-VphSynMean)/std(VphSyn(j,IdxNonZero)))<Sigma;
            if Outliar == 0
                Flag(j,i)=0;
            end
            if VrVphPth(j,i) < VrThre
                Flag(j,i)=0;
            end
        end
    end
end

end

