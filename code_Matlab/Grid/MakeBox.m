function [LonCornerTmp,LatCornerTmp,NumWgtPosTmp,NumGrNode,IdxArray]=...
    MakeBox(Lon,Lat,GrSize,WgtPosMax,IdxArray,i)

Eps=1e-4;

%% Define Coordinates for each corner
GrLonL=Lon(i,1):GrSize(i):Lon(i,2)-GrSize(i);       % Upper Left
GrLonR=Lon(i,1)+GrSize(i):GrSize(i):Lon(i,2);       % Upper Right
GrLatU=Lat(i,2):-GrSize(i):Lat(i,1)+GrSize(i);      % Lower Left
GrLatD=Lat(i,2)-GrSize(i):-GrSize(i):Lat(i,1);      % Lower Right
LonTmp=zeros(length(GrLonL),length(GrLatU),4);
LatTmp=zeros(length(GrLonL),length(GrLatU),4);
[LatTmp(:,:,1),LonTmp(:,:,1)]=meshgrid(GrLatU,GrLonL);
[LatTmp(:,:,2),LonTmp(:,:,2)]=meshgrid(GrLatU,GrLonR);
[LatTmp(:,:,3),LonTmp(:,:,3)]=meshgrid(GrLatD,GrLonL);
[LatTmp(:,:,4),LonTmp(:,:,4)]=meshgrid(GrLatD,GrLonR);
NumBox=length(GrLonL)*length(GrLatU);
LonCornerTmp=zeros(NumBox,WgtPosMax);
LatCornerTmp=zeros(NumBox,WgtPosMax);
NumWgtPosTmp=zeros(NumBox,1);
for j=1:4
    LonCornerTmp(:,j)=reshape(LonTmp(:,:,j),NumBox,1);
    LatCornerTmp(:,j)=reshape(LatTmp(:,:,j),NumBox,1);
end

%% Remove overlapping points
if     i==1
    IdxRmv=[];
elseif i==length(GrSize)
    IdxRmv=(LonCornerTmp(:,1)>Lon(i-1,1)-Eps&LonCornerTmp(:,2)<Lon(i-1,2)+Eps)&...
           (LatCornerTmp(:,1)<Lat(i-1,2)+Eps&LatCornerTmp(:,3)>Lat(i-1,1)-Eps);
end
IdxArray=max(IdxArray)+1:max(IdxArray)+NumBox-sum(IdxRmv);
LonCornerTmp(IdxRmv,:)=[];
LatCornerTmp(IdxRmv,:)=[];
NumWgtPosTmp(IdxRmv)=[];
NumGrNode=NumBox-sum(IdxRmv);

%% Additional WgtPosition
if i > 1
    for j=1:4
        switch j
            case 1; Idx1=2; Idx2=1; Idx3=4; Idx4=2; Idx5=2; % West
            case 2; Idx1=1; Idx2=2; Idx3=4; Idx4=2; Idx5=1; % East
            case 3; Idx1=4; Idx2=2; Idx3=1; Idx4=2; Idx5=4; % North
            case 4; Idx1=1; Idx2=1; Idx3=1; Idx4=2; Idx5=2; % South
        end
        if     j==1 || j==2
            Idx=(abs(LonCornerTmp(:,Idx1) - Lon(i-1,Idx2)) < Eps & ...
                    (LatCornerTmp(:,Idx3) > Lat(i-1,1) - Eps)    & ...
                    (LatCornerTmp(:,Idx4) < Lat(i-1,2) + Eps));
            LonCornerTmp(Idx,WgtPosMax)=LonCornerTmp(Idx,Idx5);
            LatCornerTmp(Idx,WgtPosMax)=(LatCornerTmp(Idx,2)+LatCornerTmp(Idx,4))/2;
        elseif j==3 || j==4
            Idx=(abs(LatCornerTmp(:,Idx1) - Lat(i-1,Idx2)) < Eps & ...
                    (LonCornerTmp(:,Idx3) > Lon(i-1,1) - Eps)    & ...
                    (LonCornerTmp(:,Idx4) < Lon(i-1,2) + Eps));
            LonCornerTmp(Idx,WgtPosMax)=(LonCornerTmp(Idx,1)+LonCornerTmp(Idx,2))/2;
            LatCornerTmp(Idx,WgtPosMax)=LatCornerTmp(Idx,Idx5);
        end
    end
end
NumWgtPosTmp(LonCornerTmp(:,WgtPosMax)~=0)=WgtPosMax;
NumWgtPosTmp(LonCornerTmp(:,WgtPosMax)==0)=WgtPosMax-1;

end
