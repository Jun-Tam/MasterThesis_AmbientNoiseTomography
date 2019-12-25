function [LonCornerTmp,LatCornerTmp,NumWgtPosTmp,NumGrNode,IdxArray]=...
    MakeBox2(Lon,Lat,GrSize,WgtPosMax,IdxArray,i)

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
    IdxRmv=(LonCornerTmp(:,1)<Lon(i,1)+GrSize(i+1)&LatCornerTmp(:,1)>Lat(i,2)-GrSize(i+1))|...
           (LonCornerTmp(:,2)>Lon(i,2)-GrSize(i+1)&LatCornerTmp(:,2)>Lat(i,2)-GrSize(i+1))|...
           (LonCornerTmp(:,3)<Lon(i,1)+GrSize(i+1)&LatCornerTmp(:,3)<Lat(i,1)+GrSize(i+1))|...
           (LonCornerTmp(:,4)>Lon(i,2)-GrSize(i+1)&LatCornerTmp(:,4)<Lat(i,1)+GrSize(i+1));
elseif i==length(GrSize)
    IdxMlt=(LonCornerTmp(:,1)>Lon(i-1,1)-Eps&LonCornerTmp(:,2)<Lon(i-1,2)+Eps)&...
           (LatCornerTmp(:,1)<Lat(i-1,2)+Eps&LatCornerTmp(:,3)>Lat(i-1,1)-Eps);
    IdxCnr=((abs(LonCornerTmp(:,1)-Lon(i-1,1))<Eps)&(abs(LatCornerTmp(:,1)-Lat(i-1,2))<Eps))|...
           ((abs(LonCornerTmp(:,2)-Lon(i-1,2))<Eps)&(abs(LatCornerTmp(:,2)-Lat(i-1,2))<Eps))|...
           ((abs(LonCornerTmp(:,3)-Lon(i-1,1))<Eps)&(abs(LatCornerTmp(:,3)-Lat(i-1,1))<Eps))|...
           ((abs(LonCornerTmp(:,4)-Lon(i-1,2))<Eps)&(abs(LatCornerTmp(:,4)-Lat(i-1,1))<Eps));       
    IdxRmv=(IdxMlt&not(IdxCnr));
else
    IdxMlt=(LonCornerTmp(:,1)>Lon(i-1,1)-Eps&LonCornerTmp(:,2)<Lon(i-1,2)+Eps)&...
           (LatCornerTmp(:,1)<Lat(i-1,2)+Eps&LatCornerTmp(:,3)>Lat(i-1,1)-Eps);
    IdxCnrIn=((abs(LonCornerTmp(:,1)-Lon(i-1,1))<Eps)&(abs(LatCornerTmp(:,1)-Lat(i-1,2))<Eps))|...
             ((abs(LonCornerTmp(:,2)-Lon(i-1,2))<Eps)&(abs(LatCornerTmp(:,2)-Lat(i-1,2))<Eps))|...
             ((abs(LonCornerTmp(:,3)-Lon(i-1,1))<Eps)&(abs(LatCornerTmp(:,3)-Lat(i-1,1))<Eps))|...
             ((abs(LonCornerTmp(:,4)-Lon(i-1,2))<Eps)&(abs(LatCornerTmp(:,4)-Lat(i-1,1))<Eps));
    IdxCnrOut=((LonCornerTmp(:,1)<Lon(i,1)+GrSize(i+1))&LatCornerTmp(:,1)>Lat(i,2)-GrSize(i+1))|...
              ((LonCornerTmp(:,2)>Lon(i,2)-GrSize(i+1))&LatCornerTmp(:,2)>Lat(i,2)-GrSize(i+1))|...
              ((LonCornerTmp(:,3)<Lon(i,1)+GrSize(i+1))&LatCornerTmp(:,3)<Lat(i,1)+GrSize(i+1))|...
              ((LonCornerTmp(:,4)>Lon(i,2)-GrSize(i+1))&LatCornerTmp(:,4)<Lat(i,1)+GrSize(i+1));         
    IdxRmv=(IdxMlt&not(IdxCnrIn)|IdxCnrOut);
end
IdxArray=max(IdxArray)+1:max(IdxArray)+NumBox-sum(IdxRmv);
LonCornerTmp(IdxRmv,:)=[];
LatCornerTmp(IdxRmv,:)=[];
NumWgtPosTmp(IdxRmv)=[];
NumGrNode=NumBox-sum(IdxRmv);

%% Additional WgtPosition
if i > 1
    for j=1:8
        switch j
            case 1; Idx1=2; Idx2=1; Idx3=4; Idx4=2; Idx5=2; % West
            case 2; Idx1=1; Idx2=2; Idx3=4; Idx4=2; Idx5=1; % East
            case 3; Idx1=4; Idx2=2; Idx3=1; Idx4=2; Idx5=4; % North
            case 4; Idx1=1; Idx2=1; Idx3=1; Idx4=2; Idx5=2; % South
            case 5; Idx1=1; Idx2=2; Idx3=3; Idx4=4; Idx5=2; % North West
            case 6; Idx1=2; Idx2=2; Idx3=4; Idx4=3; Idx5=1; % South West
            case 7; Idx1=1; Idx2=1; Idx3=1; Idx4=2; Idx5=4; % North East
            case 8; Idx1=2; Idx2=1; Idx3=2; Idx4=1; Idx5=3; % South East
        end
        if     j==1 || j==2
            Idx=(abs(LonCornerTmp(:,Idx1) - Lon(i-1,Idx2)) < Eps & ...
                    (LatCornerTmp(:,Idx3) > Lat(i-1,1) + Eps)    & ...
                    (LatCornerTmp(:,Idx4) < Lat(i-1,2) - Eps));
            LonCornerTmp(Idx,WgtPosMax-1)= LonCornerTmp(Idx,Idx5);
            LatCornerTmp(Idx,WgtPosMax-1)=(LatCornerTmp(Idx,2)+LatCornerTmp(Idx,4))/2;
        elseif j==3 || j==4
            Idx=(abs(LatCornerTmp(:,Idx1) - Lat(i-1,Idx2)) < Eps & ...
                    (LonCornerTmp(:,Idx3) > Lon(i-1,1) + Eps)    & ...
                    (LonCornerTmp(:,Idx4) < Lon(i-1,2) - Eps));
            LonCornerTmp(Idx,WgtPosMax-1)=(LonCornerTmp(Idx,1)+LonCornerTmp(Idx,2))/2;
            LatCornerTmp(Idx,WgtPosMax-1)=LatCornerTmp(Idx,Idx5);
        else
            Idx=(abs(LonCornerTmp(:,j-4) - Lon(i-1,Idx1)) < Eps & ...
                 abs(LatCornerTmp(:,j-4) - Lat(i-1,Idx2)) < Eps);
            LonCornerTmp(Idx,WgtPosMax-1)=(LonCornerTmp(Idx,Idx3)+LonCornerTmp(Idx,Idx4))/2;
            LatCornerTmp(Idx,WgtPosMax-1)= LatCornerTmp(Idx,Idx3);
            LonCornerTmp(Idx,WgtPosMax)  = LonCornerTmp(Idx,Idx5);
            LatCornerTmp(Idx,WgtPosMax)  =(LatCornerTmp(Idx,Idx4)+LatCornerTmp(Idx,Idx5))/2;
        end
    end
end
NumWgtPosTmp(LonCornerTmp(:,WgtPosMax-1)~=0&LonCornerTmp(:,WgtPosMax)~=0)=WgtPosMax;
NumWgtPosTmp(LonCornerTmp(:,WgtPosMax-1)~=0&LonCornerTmp(:,WgtPosMax)==0)=WgtPosMax-1;
NumWgtPosTmp(LonCornerTmp(:,WgtPosMax-1)==0&LonCornerTmp(:,WgtPosMax)==0)=WgtPosMax-2;

end
