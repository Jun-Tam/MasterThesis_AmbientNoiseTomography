function [VrGl,PrtbFit,Prtb0,VphSyn,Flag]=...
    GrSearch(Prtb,Freq,Vel,Dist,Amp,Data,PrtbGr,Index,VphRef0)

%% Grid search
NumGrid=length(Prtb);
VrLc=zeros(length(Prtb),1);
IdxTmp=zeros(length(Prtb),1);
for k=1:NumGrid
    WnDist=2*pi*Freq./((1+Prtb(k))*Vel')*Dist;
    Kernel=(Amp'.*besselj(0,WnDist))';
    B=pinv(Kernel'*Kernel)*Kernel'*Data;
    [~,IdxTmp(k)]=min(abs(PrtbGr-Prtb(k)));
    VrLc(k,1)=1-sum(real(B*Kernel-Data).^2)/sum(real(Data).^2);
end

if isempty(findpeaks(VrLc))==0
    [~,IdxMax]=max(VrLc);
    PrtbFit=Prtb(IdxMax);
    VphSyn=(1+PrtbFit)*mean(VphRef0(Index));
    Prtb0=PrtbFit;
    VrGl=VrLc(IdxMax);
    Flag=1;
else
    VrGl=0;
    PrtbFit=0;
    VphSyn=0;
    Prtb0=0;
    Flag=0;
end

end

