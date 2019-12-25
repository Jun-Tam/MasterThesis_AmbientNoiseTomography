function [VphMapOut,R,VphPthOut,VphErr,Fitting,Cmplxty]=...
    TmgInvCal(Grid,VphMapOut,RdErr,Weight,IdxPosWgt,IdxPosBox,...
    NumWgtPos,NumDiv,SectDist,Dist,PthDnsCvr,VphPthIn,VphPthOut,Damp)

NumNode=size(Grid,1);
NumData=length(Dist);

%% Kernel Calculation
VphErr=zeros(NumNode,1);
K=zeros(NumData,NumNode);
IdxRmv(1:NumNode)=(PthDnsCvr==0);
for i=1:NumData
    SectTime=zeros(1,NumNode);
    for j=1:NumDiv(i)
        kmax=NumWgtPos(IdxPosBox(j,i));
        IdxTmp=IdxPosWgt(j,i,1:kmax);
        SectTime(IdxTmp)=...
            SectTime(IdxTmp)+Weight(j,i,1:kmax)*SectDist(i)./VphMapOut(IdxTmp);
    end
    K(i,1:NumNode)=SectTime;
end
IdxTmp=1:NumNode;
IdxTmp(IdxRmv)=[];
K(:,IdxRmv)=[];
M=zeros(1,NumNode);
R=zeros(NumNode,NumNode);

%% Covariance Matrix
NumNodeCal=length(IdxTmp);
Cm=Damp^2*eye(NumNodeCal);

%% Inversion Process
AprErr=RdErr*Dist./(mean(VphPthIn)^2-RdErr^2);
Cd=pinv(diag(diag(AprErr*AprErr')));
dt=Dist./VphPthIn-Dist./VphPthOut;
CovPost=pinv(K'*Cd*K+Cm);
MCal=CovPost*K'*Cd*dt;
RCal=CovPost*(K'*Cd*K);
R(IdxTmp,IdxTmp)=RCal;
M(IdxTmp)=MCal;
VphErr(IdxTmp)=sqrt(diag(CovPost));
VphMapOut=VphMapOut./(1+M(1:NumNode));

%% Model velocity on used path
for i=1:NumData
    tmp=0;
    for j=1:NumDiv(i)
        kmax=NumWgtPos(IdxPosBox(j,i));
        dPath=Dist(i)/NumDiv(i);
        tmp=tmp+dPath/sum(VphMapOut(IdxPosWgt(j,i,1:kmax)).*Weight(j,i,1:kmax));
    end
    VphPthOut(i)=Dist(i)/tmp;
end
Fitting=rms(Dist./VphPthIn-Dist./VphPthOut);
Cmplxty=var(MCal);

end
