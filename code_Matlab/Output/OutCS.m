function OutCS(VsStr,VsStr0,VphRsl,Grid,GrStep,CsStep,Layer,PathOutput)

%% Make a folder
if exist(PathOutput,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOutput,'s');
end
mkdir(PathOutput);

GrLon=min(Grid(:,1)):min(GrStep):max(Grid(:,1));
GrLat=min(Grid(:,2)):min(GrStep):max(Grid(:,2));
GrLonDsp=min(Grid(:,1)):CsStep:max(Grid(:,1));
GrLatDsp=min(Grid(:,2)):CsStep:max(Grid(:,2));
[xmesh,ymesh]=meshgrid(GrLon,GrLat);
[Xmesh,Ymesh]=meshgrid(GrLonDsp,GrLatDsp);
SizeX=length(GrLonDsp);
SizeY=length(GrLatDsp);
SizeZ=length(Layer);
NumGr=SizeX*SizeY;
NumLyr=size(VsStr,2);
VsMap1=zeros(NumGr,NumLyr);
VsMap2=zeros(NumGr,NumLyr);
VsPtb1=zeros(NumGr,NumLyr);
VsPtb2=zeros(NumGr,NumLyr);
for i=1:NumLyr
    [MapTmp,~,~]=Grid2Mesh(VsStr(:,i),Grid,GrStep);
    MapTmp=interp2(xmesh,ymesh,MapTmp,Xmesh,Ymesh);
    PtbTmp=(MapTmp/VsStr0(i)-1)*100;
%     PtbTmp=(MapTmp/mean(mean(MapTmp))-1)*100;
    VsMap1(1:NumGr,i)=reshape(MapTmp',NumGr,1);
    VsMap2(1:NumGr,i)=reshape(MapTmp,NumGr,1);
    VsPtb1(1:NumGr,i)=reshape(PtbTmp',NumGr,1);
    VsPtb2(1:NumGr,i)=reshape(PtbTmp,NumGr,1);
end
VsMapTmp1=reshape(VsMap1,SizeX,SizeY,SizeZ);
VsMapTmp2=reshape(VsMap2,SizeY,SizeX,SizeZ);
VsPtbTmp1=reshape(VsPtb1,SizeX,SizeY,SizeZ);
VsPtbTmp2=reshape(VsPtb2,SizeY,SizeX,SizeZ);
[RslTmp,~,~]=Grid2Mesh(VphRsl,Grid,GrStep);
RslTmp=interp2(xmesh,ymesh,RslTmp,Xmesh,Ymesh);
VphRslTmp1=reshape(RslTmp',SizeX,SizeY);
VphRslTmp2=reshape(RslTmp,SizeY,SizeX);

%% Elevation Input
Fid=fopen('elv.txt','r');
Data=textscan(Fid,'%f %f %f\n');
ElvLon=Data{1};
ElvLat=Data{2};
Elevation=Data{3};
[xq,yq]=meshgrid(GrLonDsp,GrLatDsp);
F=TriScatteredInterp(ElvLon,ElvLat,Elevation,'natural'); 
AltMap=F(xq,yq);

%% Post Shock
% 3 month after the main shock (Determined by Okada et al., 2012)
% Fid=fopen('hypo.txt','r');
% Dataset=textscan(Fid,repmat('%f ',1,13),'HeaderLines',1);
% fclose(Fid);
% LatPS=Dataset{7};
% LonPS=Dataset{8};
% DepPS=Dataset{9};
% MagPS=Dataset{10};

% 2008-2012 (Determined by Okada et al., 2012)
Fid=fopen('hypo.all.txt','r');
Dataset=textscan(Fid,['%s',repmat('%f ',1,14)]);
fclose(Fid);
LatPS1=Dataset{8};
LonPS1=Dataset{9};
DepPS1=Dataset{10};
MagPS1=Dataset{11};

% 2008-2012 (Determined by Ambiet noise tomography)
% Fid=fopen('tomoDD.reloc.2008-2010.txt','r');
% Dataset=textscan(Fid,repmat('%f ',1,24));
% fclose(Fid);
% LatPS1=Dataset{2};
% LonPS1=Dataset{3};
% DepPS1=Dataset{4};
% MagPS1=Dataset{17};

% 2008-2012 (Determined by Ambiet noise tomography)
Fid=fopen('tomoDD.reloc.2011-2012.txt','r');
Dataset=textscan(Fid,repmat('%f ',1,24));
fclose(Fid);
LatPS2=Dataset{2};
LonPS2=Dataset{3};
DepPS2=Dataset{4};
MagPS2=Dataset{17};

%% Volcano position
Fid=fopen('VlcQ.txt','r');
VlcInfo=textscan(Fid,'%f %f\n');
VlcLon=VlcInfo{1};
VlcLat=VlcInfo{2};
fclose(Fid);

%% Reflector
Fid=fopen('RefInfo.txt','r');
DataSet=cell2mat(textscan(Fid,repmat('%f ',1,8)));
fclose(Fid);
RefLon=DataSet(:,1);
RefLat=DataSet(:,3);
RefDep=DataSet(:,5);
RefStr=DataSet(:,7);
RefDip=DataSet(:,8);

%% Station
% Fid=fopen('StationPresent.txt','r'); % ANT
% DataSet=cell2mat(textscan(Fid,repmat('%f ',1,2)));
% fclose(Fid);
% StLon=DataSet(:,1);
% StLat=DataSet(:,2);

Fid=fopen('Station_Okada2012.txt','r'); % Okada et al. (2012)
DataSet=textscan(Fid,'%s %f %f'); 
fclose(Fid);
StLon=cell2mat(DataSet(:,3));
StLat=cell2mat(DataSet(:,2));

%% Cross-section
OutCS2(VlcLat,VlcLon,'NS/',GrLonDsp,GrLatDsp,Grid(:,1),Grid(:,2),Layer,LonPS1,...
    LatPS1,DepPS1,MagPS1,LonPS2,LatPS2,DepPS2,MagPS2,AltMap,VsMapTmp1,VsPtbTmp1,...
    VphRslTmp1,RefLat,RefLon,RefDep,RefStr,RefDip,StLon,StLat,PathOutput);
OutCS2(VlcLon,VlcLat,'EW/',GrLatDsp,GrLonDsp,Grid(:,2),Grid(:,1),Layer,LatPS1,...
    LonPS1,DepPS1,MagPS1,LatPS2,LonPS2,DepPS2,MagPS2,AltMap',VsMapTmp2,VsPtbTmp2,...
    VphRslTmp2,RefLon,RefLat,RefDep,RefStr,RefDip,StLat,StLon,PathOutput);

end

