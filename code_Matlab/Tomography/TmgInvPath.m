function [Weight,IdxPosWgt,IdxPosBox,SplitNum,SectDist,PthDnsCvr,PthAzmCvr]=...
    TmgInvPath(StPos,Azim,Dist,IdxCorner,GrStep,Grid,NumWgtPos,NumSect,...
    NumBinAzm,FlagSect)

%% Initial Setting
DataLength=length(Dist);
NumNode=size(Grid,1);
DelLon=zeros(1,DataLength);
DelLat=zeros(1,DataLength);
SectDist=zeros(1,DataLength);
SplitNum=zeros(1,DataLength);
PathPtLon=zeros(DataLength,NumSect);
PathPtLat=zeros(DataLength,NumSect);
Weight=zeros(NumSect,DataLength,max(NumWgtPos));
IdxPosWgt=zeros(NumSect,DataLength);
IdxPosBox=zeros(NumSect,DataLength);
PthDnsCvr=zeros(1,NumNode);
PthAzmCvr=zeros(1,NumNode);
AzmCnt=zeros(NumBinAzm,NumNode);
BinSize=180/(NumBinAzm);
Bin=zeros(NumBinAzm,2);
Bin(1:NumBinAzm-1,1)=BinSize/2:BinSize:180-BinSize;
Bin(1:NumBinAzm-1,2)=BinSize+BinSize/2:BinSize:180-BinSize/2;

%% Path Processing
for i=1:DataLength
    if FlagSect==0
        SplitNum(i)=NumSect;
    elseif FlagSect==1
        SplitNum(i)=round(Dist(i))*NumSect;
    end
    SectDist(i)=Dist(i)/SplitNum(i);
    DelLon(i)=(StPos(i,3)-StPos(i,1))/(2*SplitNum(i));
    DelLat(i)=(StPos(i,4)-StPos(i,2))/(2*SplitNum(i));
    for j=1:SplitNum(i)
        PathPtLon(i,j)=StPos(i,1)+(2*j-1)*DelLon(i);
        PathPtLat(i,j)=StPos(i,2)+(2*j-1)*DelLat(i);
        IdxPosBox(j,i)= find(Grid(IdxCorner(:,1),1) <= PathPtLon(i,j) & ...
                             Grid(IdxCorner(:,4),1) >  PathPtLon(i,j) & ...
                             Grid(IdxCorner(:,1),2) >= PathPtLat(i,j) & ...
                             Grid(IdxCorner(:,4),2) <  PathPtLat(i,j));     
        IdxTmp=1:NumWgtPos(IdxPosBox(j,i));
        IdxPosWgt(j,i,IdxTmp)=IdxCorner(IdxPosBox(j,i),IdxTmp);
        DstLonCnr=(Grid(IdxCorner(IdxPosBox(j,i),IdxTmp),1)-PathPtLon(i,j)).^2;
        DstLatCnr=(Grid(IdxCorner(IdxPosBox(j,i),IdxTmp),2)-PathPtLat(i,j)).^2;        
        DstCnr=sqrt(DstLonCnr+DstLatCnr)/GrStep(IdxPosBox(j,i));
        WeightTmp=interp1([0 sqrt(2)/2 sqrt(2)],[1 0.25 0],DstCnr);
        Weight(j,i,IdxTmp)=WeightTmp/sum(WeightTmp);
        PthDnsCvr(IdxCorner(IdxPosBox(j,i),IdxTmp))=...
            PthDnsCvr(IdxCorner(IdxPosBox(j,i),IdxTmp))+WeightTmp'/sum(WeightTmp);
        IdxBin=find(Azim(i) > Bin(:,1) & Azim(i) < Bin(:,2));
        if isempty(IdxBin)==1
            IdxBin=NumBinAzm;
        end
        AzmCnt(IdxBin,IdxCorner(IdxPosBox(j,i),IdxTmp))=...
            AzmCnt(IdxBin,IdxCorner(IdxPosBox(j,i),IdxTmp))+1;
    end
end
PthDnsCvr=PthDnsCvr*DataLength/sum(round(Dist)*NumSect);
for i=1:NumNode
    if sum(AzmCnt(:,i)) > 100
        Numer=sum(AzmCnt(:,i)./sum(AzmCnt(:,i)));
        Denom=NumBinAzm*(max(AzmCnt(:,i))/sum(AzmCnt(:,i)));
        PthAzmCvr(i)=Numer/Denom;
    end
end

end
