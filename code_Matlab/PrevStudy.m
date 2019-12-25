%% Initial Setting
clear all;
Path='/Users/Jun/Dropbox';
PathHome=strcat(Path,'/Tomography');
PathIn=strcat(PathHome,'/Input/');
PathOut=strcat(PathHome,'/Output/');
addpath(genpath(pwd));
addpath(genpath(strcat(Path,'/Code/Matlab')));
addpath(genpath(strcat(Path,'/InfoSource')));
cd(PathHome);

%% Data input
Fid=fopen('grid.result.txt','r');
Dataset=cell2mat(textscan(Fid,repmat('%f ',1,9),'HeaderLines',1));
Layer=[ones(8,1);2;];
Depth=[0;cumsum(Layer(1:end-1))];
GrSize=0.04;
LatRange=[38.54, 39.10];
LonRange=[140.36,141.08];
[~,Grid,~,~]=RglGrid(LonRange,LatRange,GrSize);
GrLon=min(Grid(:,1)):0.04:max(Grid(:,1));
GrLat=min(Grid(:,2)):0.04:max(Grid(:,2));
NumNode=length(GrLon)*length(GrLat);
Dep=unique(Dataset(:,3));
NumDepIn=length(Dep);
NumDep=length(Depth);
NumGrPt=length(Dataset(:,1));
Lon=reshape(Dataset(:,1),NumGrPt/NumDepIn,NumDepIn);
Lat=reshape(Dataset(:,2),NumGrPt/NumDepIn,NumDepIn);
Vp=reshape(Dataset(:,4),NumGrPt/NumDepIn,NumDepIn);
Vs=reshape(Dataset(:,5),NumGrPt/NumDepIn,NumDepIn);
VpVs=reshape(Dataset(:,6),NumGrPt/NumDepIn,NumDepIn);
Res=reshape(Dataset(:,9),NumGrPt/NumDepIn,NumDepIn);

%% JMA Model
Fid=fopen('vjma2001','r');
Dataset=textscan(Fid,repmat('%f ',1,3));
fclose(Fid);
VpJMA=interp1(Dataset{3},Dataset{1},Depth);
VsJMA=interp1(Dataset{3},Dataset{2},Depth);

%% Interpolation
[xq,yq]=meshgrid(GrLon,GrLat);
VpMapTmp=zeros(NumNode,1);
VsMapTmp=zeros(NumNode,1);
VpVsMapTmp=zeros(NumNode,1);
ResMapTmp=zeros(NumNode,1);
VpMap=zeros(NumNode,NumDep);
VsMap=zeros(NumNode,NumDep);
VpVsMap=zeros(NumNode,NumDep);
ResMap=zeros(NumNode,NumDep);
for i=1:length(Dep)
    FVp=TriScatteredInterp(Lon(:,i),Lat(:,i),Vp(:,i),'natural');
    FVs=TriScatteredInterp(Lon(:,i),Lat(:,i),Vs(:,i),'natural'); 
    FVpVs=TriScatteredInterp(Lon(:,i),Lat(:,i),VpVs(:,i),'natural'); 
    FRes=TriScatteredInterp(Lon(:,i),Lat(:,i),Res(:,i),'natural');
    VpMapTmp(:,i)=reshape(FVp(xq,yq),NumNode,1);
    VsMapTmp(:,i)=reshape(FVs(xq,yq),NumNode,1);
    VpVsMapTmp(:,i)=reshape(FVpVs(xq,yq),NumNode,1);
    ResMapTmp(:,i)=reshape(FRes(xq,yq),NumNode,1);
end
for i=1:NumNode
    VpMap(i,:)=interp1(Dep,VpMapTmp(i,:),Depth,'linear','extrap');
    VsMap(i,:)=interp1(Dep,VsMapTmp(i,:),Depth,'linear','extrap');
    VpVsMap(i,:)=interp1(Dep,VpVsMapTmp(i,:),Depth,'linear','extrap');
    ResMap(i,:)=interp1(Dep,ResMapTmp(i,:),Depth,'linear','extrap');
end

%% Output
PathOutput=strcat(PathOut,'Ref/');
if exist(PathOutput,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOutput,'s');
end
mkdir(PathOutput);

for i=1:NumDep
    Fid=fopen(strcat(PathOutput,num2str(Depth(i),'%0.1f'),'km.txt'),'w');
    IdxNonNan=isnan(VpMap(:,i))~=1;
    VpPtb=zeros(NumNode,1);
    VsPtb=zeros(NumNode,1);
    VpPtb(IdxNonNan)=(VpMap(IdxNonNan,i)/VpJMA(i)-1)*100;
    VsPtb(IdxNonNan)=(VsMap(IdxNonNan,i)/VsJMA(i)-1)*100;
    for j=1:NumNode
        fprintf(Fid(1),[repmat('%f ',1,8),'\n'],Grid(j,1:2),VpMap(j,i),...
            VpPtb(j),VsMap(j,i),VsPtb(j),VpVsMap(j,i),ResMap(j,i));
    end
    fclose(Fid);
end


