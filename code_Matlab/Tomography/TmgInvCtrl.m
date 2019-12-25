function [VphPthIn,VphPthOut,VphMapIn,VphMapOut,VphErr,VphRsl,VphCmp,VphFit,...
    PthDnsCvr,PthAzmCvr]=TmgInvCtrl(RdErr,Grid,GrStep,IdxCorner,NumWgtPos,...
    Damp,FlagMdl,NumSect,NumBinAzm,FlagSect,VphPthIn,StAzim,StDist,StPos,VphRef)

%% Path Processing
[Weight,IdxPosWgt,IdxPosBox,NumDiv,SectDist,PthDnsCvr,PthAzmCvr]=...
    TmgInvPath(StPos,StAzim,StDist,IdxCorner,GrStep,Grid,NumWgtPos,NumSect,...
    NumBinAzm,FlagSect);

%% Initial Model
[VphPthIn,VphMap,VphPth,VphMapIn]=TmgInvPrep(VphPthIn,VphRef,StDist,Grid,...
    GrStep,NumWgtPos,Weight,IdxPosWgt,IdxPosBox,NumDiv,FlagMdl);

%% Inversion
NumDamp=length(Damp);
NumPath=length(VphPth);
NumNode=length(VphMap);
VphPthOut=zeros(NumPath,NumDamp);
VphMapOut=zeros(NumNode,NumDamp);
VphErr=zeros(NumNode,NumDamp);
VphRsl=zeros(NumNode^2,NumDamp);
VphCmp=zeros(NumDamp,1);
VphFit=zeros(NumDamp,1);
for i=1:NumDamp
    [VphMapOut(:,i),R,VphPthOut(:,i),VphErr(:,i),VphFit(i),VphCmp(i)]=...
        TmgInvCal(Grid,VphMap,RdErr,Weight,IdxPosWgt,IdxPosBox,NumWgtPos,...
        NumDiv,SectDist,StDist,PthDnsCvr,VphPthIn,VphPth,Damp(i));
    VphRsl(:,i)=reshape(R,NumNode^2,1);
end

end
