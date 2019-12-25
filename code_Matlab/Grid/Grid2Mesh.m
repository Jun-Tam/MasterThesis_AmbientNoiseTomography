function [MeshImage,GrLon,GrLat]=Grid2Mesh(DataGrid,Grid,GrStep)

GrLon=min(Grid(:,1)):min(GrStep):max(Grid(:,1));
GrLat=min(Grid(:,2)):min(GrStep):max(Grid(:,2));
[LonMesh,LatMesh]=meshgrid(GrLon,GrLat);
F=TriScatteredInterp(Grid(:,1),Grid(:,2),DataGrid,'natural');
MeshImage=F(LonMesh,LatMesh);

end