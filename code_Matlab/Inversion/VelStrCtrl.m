function [VsStr,VsErr,VsRsl,VphStr,VphErr,VsCmp,VsFit,IdxItr,RmsItr,IdxDmp,Flag]=...
    VelStrCtrl(VsStrIn,VphMap,VphErrIn,DmpStr,ItrStr,CrDist,Layer,Freq,PathK,IdxDmpMdl)

%% Initial setting
if isempty(IdxDmpMdl)==1
    NumDmp=length(DmpStr);
else
    NumDmp=1;
end
NumData=size(VphMap,2);
NumLyr=length(Layer);
NumFreq=length(Freq);
VsStr=zeros(NumData,NumLyr);
VsErr=zeros(NumData,NumLyr);
VsRsl=zeros(NumData,NumLyr);
VphStr=zeros(NumData,NumFreq);
VphErr=zeros(NumData,NumFreq);
VsFit=zeros(NumData,NumDmp);
VsCmp=zeros(NumData,NumDmp);
IdxItr=zeros(1,NumData);
IdxDmp=zeros(1,NumData);
RmsItrTmp=zeros(NumData,NumDmp,ItrStr);
RmsItr=zeros(NumData,ItrStr);
Flag=false(NumData,NumDmp);

%% Calculation
for i=1:NumData
    VsStrTmp=zeros(NumLyr,NumDmp);
    VsErrTmp=zeros(NumLyr,NumDmp);
    VsRslTmp=zeros(NumLyr,NumDmp);
    VphStrTmp=zeros(NumFreq,NumDmp);
    VphErrTmp=zeros(NumFreq,NumDmp);
    IdxItrTmp=zeros(1,NumDmp);
    if isempty(IdxDmpMdl)==0
        j=1;
        [VsStrTmp(:,j),VsErrTmp(:,j),VsRslTmp(:,j),VsFit(i,j),VsCmp(i,j),...
            VphStrTmp(:,j),VphErrTmp(:,j),RmsItrTmp(i,j,:),IdxItrTmp(j)]=...
            VelStrCal(VsStrIn,VphMap(:,i),VphErrIn,DmpStr(IdxDmpMdl(i)),...
            ItrStr,CrDist,Layer,Freq,PathK);        
    else
        for j=1:NumDmp
            [VsStrTmp(:,j),VsErrTmp(:,j),VsRslTmp(:,j),VsFit(i,j),VsCmp(i,j),...
                VphStrTmp(:,j),VphErrTmp(:,j),RmsItrTmp(i,j,:),IdxItrTmp(j)]=...
                VelStrCal(VsStrIn,VphMap(:,i),VphErrIn,DmpStr(j),ItrStr,...
                CrDist,Layer,Freq,PathK);
        end
    end
    Min=VsCmp(i,NumDmp);
    Max=VsFit(i,NumDmp);
    Flag(i,NumDmp)=true;
    for j=NumDmp-1:-1:1
        if VsFit(i,j) < Max && VsCmp(i,j) > Min
            Max=VsFit(i,j);
            Min=VsCmp(i,j);
            Flag(i,j)=true;
        end
    end
    PtOpt=FindOptPt(VsCmp(i,Flag(i,:)),VsFit(i,Flag(i,:)),DmpStr(Flag(i,:)));
    IdxDmp(i)=max(find(Flag(i,:)==1,PtOpt));
    VsStr(i,:)=VsStrTmp(:,IdxDmp(i));
    VsErr(i,:)=VsErrTmp(:,IdxDmp(i));
    VsRsl(i,:)=VsRslTmp(:,IdxDmp(i));
    VphStr(i,:)=VphStrTmp(:,IdxDmp(i));
    VphErr(i,:)=VphErrTmp(:,IdxDmp(i));
    RmsItr(i,:)=RmsItrTmp(i,IdxDmp(i),:);
    IdxItr(i)=IdxItrTmp(IdxDmp(i));
end

end
