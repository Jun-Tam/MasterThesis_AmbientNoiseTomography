function [VphRef0,VrVphRef0,VrRef0,Amp,FreqRef,SpRef,SpReal,SpImag]=...
    RefVelCal(WfsFreq,Dist,VelGrRange,VelGrStep,FreqRange,FreqStep,Fs,...
    WfLength,PrtbRef)

%% Setting
Freq=Fs/2*linspace(-1,1,WfLength);
VelGr=min(VelGrRange):VelGrStep:max(VelGrRange);
IdxSpP=find(Freq > min(FreqRange) - 2*FreqStep & ...
            Freq < max(FreqRange) + 2*FreqStep);
IdxSpN=sort(find(Freq < -min(FreqRange) + 2*FreqStep & ...
                 Freq > -max(FreqRange) - 2*FreqStep),'descend');
NumFreq=length(IdxSpP);
NumGrid=length(VelGr);
NumData=length(Dist);
A=zeros(NumFreq,NumGrid);
Amp=zeros(NumFreq,1);
SpReal=zeros(2*NumData,NumFreq);
SpImag=zeros(2*NumData,NumFreq);
VrRef0=zeros(NumFreq,NumGrid);
VphRef0=zeros(NumFreq,1);
VrVphRef0=zeros(NumFreq,1);
FreqRef=Freq(IdxSpP);
SpRef=WfsFreq(:,IdxSpP);

%% Grid Search
for i=1:NumFreq
    for j=1:NumGrid
        WaveNum=(2*pi)*Freq(IdxSpP(i))/VelGr(j);
        Kernel=besselj(0,WaveNum*Dist);
        A(i,j)=Kernel\real(SpRef(:,i));
        Misfit=Kernel*reshape(A(i,j),1,1) - real(SpRef(:,i));
        VrRef0(i,j)=1-sum(abs(Misfit).^2)/sum(real(SpRef(:,i)).^2);
    end
    if i == 1
        IdxRange=1:length(NumGrid);
        IdxBase=0;
    else
        IdxRange=VelGr > VphRef0(i-1)*(1 - PrtbRef) & ...
                 VelGr < VphRef0(i-1)*(1 + PrtbRef);
        IdxBase=find(IdxRange, 1 )-1;        
    end
    [~,IdxMax]=max(VrRef0(i,IdxRange));
    Amp(i)=A(i,IdxMax+IdxBase);
    VphRef0(i)=VelGr(IdxMax+IdxBase);
    VrVphRef0(i)=max(VrRef0(i,IdxRange));
end

%% For Display
for i=1:NumData
    SpReal(i,1:NumFreq)=real(WfsFreq(i,IdxSpN));
    SpImag(i,1:NumFreq)=imag(WfsFreq(i,IdxSpN));
    SpReal(i+NumData,1:NumFreq)=real(WfsFreq(i,IdxSpP));
    SpImag(i+NumData,1:NumFreq)=imag(WfsFreq(i,IdxSpP));
end

end
