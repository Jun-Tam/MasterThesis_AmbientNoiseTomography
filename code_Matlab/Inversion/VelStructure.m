function [VphStr,VphErr,VsStr,VsErr,VsRsl,VsCmp,VsFit,IdxDmp,IdxItr,RmsItr,...
    VsStr0,FlagDmp]=VelStructure(VphMap,VphErrIn,PthDnsCvr,Layer,Grid,GrStep,...
    Freq,ItrStr,CrDist,DmpStr,VphRef,FlagVs0,PathK,IdxDmpMdl,FlagRcv)

%% Initial Setting
if isempty(IdxDmpMdl)==0 && FlagRcv==2
    NumDmp=1;
else
    NumDmp=length(DmpStr);    
end

NumNode=size(Grid,1);
NumFreq=length(Freq);
NumLyr=length(Layer);
VphStr=zeros(NumNode,NumFreq);
VphErr=zeros(NumNode,NumFreq);
VsStr=zeros(NumNode,NumLyr);
VsErr=zeros(NumNode,NumLyr);
VsRsl=zeros(NumNode,NumLyr);
VsCmp=zeros(NumNode,NumDmp);
VsFit=zeros(NumNode,NumDmp);
IdxDmp=zeros(1,NumNode);
FlagDmp=zeros(NumNode,NumDmp);
IdxItr=zeros(1,NumNode);
RmsItr=zeros(NumNode,ItrStr);
MinPthDns=min(PthDnsCvr,[],2);
IdxPath=find(MinPthDns > 0);
IdxNoPath=find(MinPthDns == 0);

%% Initial Model
VsStr0=VelStrPrep(Layer,Freq,VphRef,VphErrIn,CrDist,DmpStr,ItrStr,FlagVs0,PathK);

tic;
%% Velocity Structure (Path Density < 1)
if isempty(IdxDmpMdl)==0 && FlagRcv==2
    DmpStrTmp=DmpStr(mode(IdxDmpMdl(IdxNoPath)));
    [VsStrOutTmp,VsErrTmp,VsRslTmp,VphStrTmp,VphErrTmp,VsCmpTmp,VsFitTmp,...
        IdxItrTmp,RmsItrTmp,DmpIdxTmp,FlagDmpTmp]=VelStrCtrl(VsStr0',...
        mean(VphMap(IdxNoPath,:),1)',VphErrIn,DmpStrTmp,ItrStr,CrDist,...
        Layer,Freq,PathK,[]);
else
    [VsStrOutTmp,VsErrTmp,VsRslTmp,VphStrTmp,VphErrTmp,VsCmpTmp,VsFitTmp,...
        IdxItrTmp,RmsItrTmp,DmpIdxTmp,FlagDmpTmp]=VelStrCtrl(VsStr0',...
        mean(VphMap(IdxNoPath,:),1)',VphErrIn,DmpStr,ItrStr,CrDist,...
        Layer,Freq,PathK,[]);    
end
for i=1:length(IdxNoPath)
    j=IdxNoPath(i);
    VsStr(j,:)=VsStrOutTmp;
    VsRsl(j,:)=VsRslTmp;
    VsErr(j,:)=VsErrTmp;
    VphStr(j,:)=VphStrTmp;
    VphErr(j,:)=VphErrTmp;
    VsFit(j,:)=VsFitTmp;
    VsCmp(j,:)=VsCmpTmp;
    IdxItr(j)=IdxItrTmp;
    RmsItr(j,:)=RmsItrTmp;
    IdxDmp(j)=DmpIdxTmp;
    FlagDmp(j,:)=FlagDmpTmp;
end

%% Velocity Structure (Path Density >= 1)
if isempty(IdxDmpMdl)==0 && FlagRcv==2
    [VsStr(IdxPath,:),VsErr(IdxPath,:),VsRsl(IdxPath,:),VphStr(IdxPath,:),...
        VphErr(IdxPath,:),VsCmp(IdxPath,:),VsFit(IdxPath,:),IdxItr(IdxPath),...
        RmsItr(IdxPath,:),IdxDmp(IdxPath),FlagDmp(IdxPath,:)]=...
        VelStrCtrl(VsStr0',VphMap(IdxPath,:)',VphErrIn,...
        DmpStr,ItrStr,CrDist,Layer,Freq,PathK,IdxDmpMdl(IdxPath));
else    
    [VsStr(IdxPath,:),VsErr(IdxPath,:),VsRsl(IdxPath,:),VphStr(IdxPath,:),...
        VphErr(IdxPath,:),VsCmp(IdxPath,:),VsFit(IdxPath,:),IdxItr(IdxPath),...
        RmsItr(IdxPath,:),IdxDmp(IdxPath),FlagDmp(IdxPath,:)]=...
        VelStrCtrl(VsStr0',VphMap(IdxPath,:)',VphErrIn,...
        DmpStr,ItrStr,CrDist,Layer,Freq,PathK,[]);
end
LapseTime=toc;
display(strcat('Velocity Structure:',num2str(LapseTime),'sec'));

end
