function DataGrid=Mesh2Grid(DataMesh,Grid,GrStep)

GrLon=min(Grid(:,1)):min(GrStep):max(Grid(:,1));
GrLat=min(Grid(:,2)):min(GrStep):max(Grid(:,2));
[LonMesh,LatMesh]=meshgrid(GrLon,GrLat);
DataGrid=interp2(LonMesh,LatMesh,DataMesh,Grid(:,1),Grid(:,2));

end

