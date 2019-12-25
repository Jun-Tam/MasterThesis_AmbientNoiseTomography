function [GrStep,Grid,IdxCorner,NumWgtPos] = RglGrid(Lon,Lat,GrSize)

%% Grid make
[GrLon,GrLat]=meshgrid(Lon(1):GrSize:Lon(2),Lat(2):-GrSize:Lat(1));
NumGrPt=numel(GrLon);
Grid=zeros(NumGrPt,2);
Grid(:,1)=reshape(GrLon',NumGrPt,1);
Grid(:,2)=reshape(GrLat',NumGrPt,1);

%% Box assign
NumLon=size(GrLon,2);
NumLat=size(GrLon,1);
NumBox=(NumLon-1)*(NumLat-1);
NumWgtPos=repmat(4,NumBox,1);
GrStep=repmat(GrSize,NumBox,1);
IdxCorner=zeros(NumBox,4);
for i=1:NumLat-1
    IdxArray=(i-1)*(NumLon-1)+1:i*(NumLon-1);
    IdxCorner(IdxArray,1)=(i-1)*NumLon+1:i*NumLon-1;
    IdxCorner(IdxArray,2)=(i-1)*NumLon+2:i*NumLon;
    IdxCorner(IdxArray,3)=i*NumLon+1:(i+1)*NumLon-1;
    IdxCorner(IdxArray,4)=i*NumLon+2:(i+1)*NumLon;
end

end

