function DataCmp= NanCmp(Data,Grid,GrStep,IdxNaN)

[DataMesh,~,~]=Grid2Mesh(Data,Grid,GrStep);
DataGrid=Mesh2Grid(inpaint_nans(DataMesh),Grid,GrStep);
DataCmp=DataGrid(IdxNaN);

end

