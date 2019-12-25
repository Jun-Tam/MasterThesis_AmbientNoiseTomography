function [VphPthIn,VphPthOut,VphMapIn,VphMapOut,VphRms,VphErr,VphRsl,...
    VphCmp,VphFit,DmpOpt,PthDnsCvr,PthAzmCvr,Grid,GrStep,Cmplxty,Fitting]=...
    Tomography(VphRef,VphPthMdlIn,StAzim,StDist,StPos,Freq,DmpTmg,RdErr,...
    FlagMdl,Lon,Lat,GrSize,NumSect,NumBinAzm,FlagSect,NumPathMax)

tic;
%% Mesh grid definition
[GrStep,Grid,IdxCorner,NumWgtPos]=RglGrid(Lon,Lat,GrSize);

%% Initial Setting
NumNode=size(Grid,1);
NumFreq=length(Freq);
NumDamp=length(DmpTmg);
VphMapIn=zeros(NumNode,NumFreq);
VphMapOut=zeros(NumNode,NumFreq);
VphRsl=zeros(NumNode,NumFreq);
VphRms=zeros(1,NumFreq);
VphCmp=zeros(1,NumFreq);
VphFit=zeros(1,NumFreq);
DmpOpt=zeros(1,NumFreq);
Cmplxty=zeros(NumDamp,NumFreq);
Fitting=zeros(NumDamp,NumFreq);
VphPthIn=zeros(NumPathMax,NumFreq);
VphPthOut=zeros(NumPathMax,NumFreq);
VphErr=zeros(NumNode,1);
PthDnsCvr=zeros(NumNode,NumFreq);
PthAzmCvr=zeros(NumNode,NumFreq);
ApriErr=mean(RdErr);

%% Find Optimal Damping parameter
for j=1:NumFreq
    IdxNonZero=VphPthMdlIn(:,j)~=0;
    VphPthIn0=VphPthMdlIn(IdxNonZero,j);
    Dist0=StDist(IdxNonZero,j);
    Azim0=StAzim(IdxNonZero,j);
    Pos0=StPos(IdxNonZero,1:4,j);
    [~,~,~,~,~,~,Cmplxty(:,j),Fitting(:,j),~,~]=...
        TmgInvCtrl(ApriErr,Grid,GrStep,IdxCorner,NumWgtPos,DmpTmg,FlagMdl,...
        NumSect,NumBinAzm,FlagSect,VphPthIn0,Azim0,Dist0,Pos0,VphRef(j));
    DmpOpt(j)=DmpTmg(FindOptPt(Cmplxty(:,j),Fitting(:,j),DmpTmg));
end

%% Tomographic inversion
for j=1:NumFreq
    IdxNonZero=VphPthMdlIn(:,j)~=0;
    VphPthIn0=VphPthMdlIn(IdxNonZero,j);
    Dist0=StDist(IdxNonZero,j);
    Azim0=StAzim(IdxNonZero,j);
    Pos0=StPos(IdxNonZero,1:4,j);
    [VphPthInTmp,VphPthOutTmp,VphMapIn(:,j),VphMapOut(:,j),VphErr(:,j),R,...
        VphCmp(j),VphFit(j),PthDnsCvr(:,j),PthAzmCvr(:,j)]=...
        TmgInvCtrl(ApriErr,Grid,GrStep,IdxCorner,NumWgtPos,DmpOpt(j),FlagMdl,...
        NumSect,NumBinAzm,FlagSect,VphPthIn0,Dist0,Azim0,Pos0,VphRef(j));
    NumPath=length(VphPthInTmp);
    VphPthIn(1:NumPath,j)=VphPthInTmp;
    VphPthOut(1:NumPath,j)=VphPthOutTmp;
    VphRms(j)=rms(VphPthIn(:,j)-VphPthOut(:,j));
    VphRsl(:,j)=diag(reshape(R,NumNode,NumNode));
end
LapseTime=toc;
display(strcat('Tomography:',num2str(LapseTime),'sec'));

end
