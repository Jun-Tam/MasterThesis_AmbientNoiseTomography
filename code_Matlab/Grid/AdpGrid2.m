function [GrStep,Grid,IdxCorner,NumWgtPos] = AdpGrid2(Lon,Lat,GrSize)


%% Make Grid
Grid=MakeGrid2(Lon,Lat,GrSize);

%% Initial Setting
if length(GrSize)==1
    WgtPosMax=4;
    GrLonL=Lon(1):+GrSize:Lon(2)-GrSize;       % Upper Left
    GrLonR=Lon(1)+GrSize:+GrSize:Lon(2);       % Upper Right
    GrLatU=Lat(2):-GrSize:Lat(1)+GrSize;       % Lower Left
    GrLatD=Lat(2)-GrSize:-GrSize:Lat(1);       % Lower Right
    LonTmp=zeros(length(GrLonL),length(GrLatU),4);
    LatTmp=zeros(length(GrLonL),length(GrLatU),4);
    [LatTmp(:,:,1),LonTmp(:,:,1)]=meshgrid(GrLatU,GrLonL);
    [LatTmp(:,:,2),LonTmp(:,:,2)]=meshgrid(GrLatU,GrLonR);
    [LatTmp(:,:,3),LonTmp(:,:,3)]=meshgrid(GrLatD,GrLonL);
    [LatTmp(:,:,4),LonTmp(:,:,4)]=meshgrid(GrLatD,GrLonR);
    NumBox=length(GrLonL)*length(GrLatU);
    LonCornerTmp=zeros(NumBox,WgtPosMax);
    LatCornerTmp=zeros(NumBox,WgtPosMax);
    IdxCorner=zeros(NumBox,WgtPosMax);    
    for j=1:WgtPosMax
        LonCornerTmp(:,j)=reshape(LonTmp(:,:,j),NumBox,1);
        LatCornerTmp(:,j)=reshape(LatTmp(:,:,j),NumBox,1);
    end    
    for j=1:NumBox
        for k=1:WgtPosMax
            Residual=abs(Grid(:,1)-LonCornerTmp(j,k))+...
                     abs(Grid(:,2)-LatCornerTmp(j,k));
            IdxCorner(j,k)=find(Residual==min(Residual));
        end
    end
    NumWgtPos=WgtPosMax*ones(NumBox,1);
    GrStep=GrSize*ones(NumBox,1);
else
    WgtPosMax=6;
    IdxArray=0;
    IdxCorner=zeros(1000,WgtPosMax);
    NumWgtPos=zeros(1000,1);
    GrStep=zeros(1000,1);
    
    %% Assign Index to 4 Corners
    for i=1:length(GrSize)
        [LonCornerTmp,LatCornerTmp,NumWgtPosTmp,NumGrNode,IdxArray]=...
            MakeBox2(Lon,Lat,GrSize,WgtPosMax,IdxArray,i);
        IdxCornerTmp=zeros(NumGrNode,WgtPosMax);
        for j=1:NumGrNode
            for k=1:NumWgtPosTmp(j)
                Residual=abs(Grid(:,1)-LonCornerTmp(j,k))+...
                         abs(Grid(:,2)-LatCornerTmp(j,k));
                IdxCornerTmp(j,k)=find(Residual==min(Residual));
            end
        end
        IdxCorner(IdxArray,:)=IdxCornerTmp;
        NumWgtPos(IdxArray,1)=NumWgtPosTmp;
        GrStep(IdxArray)=GrSize(i);
    end
    IdxRmv=(GrStep==0);
    IdxCorner(IdxRmv,:)=[];
    NumWgtPos(IdxRmv)=[];
    GrStep(IdxRmv)=[];    
end

end

