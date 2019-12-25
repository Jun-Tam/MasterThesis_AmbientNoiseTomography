%% Body-wave seismic tomography
Fid=fopen('vel_Okada2012.txt','r');
Dataset=cell2mat(textscan(Fid,repmat('%f ',1,7),'HeaderLines',1));
Dep=unique(Dataset(:,3));
Depth=[0;cumsum(Layer(1:end-1))];
NumDep=length(Dep);
NumGrPt=length(Dataset(:,1));
% NumNode=NumGrPt/NumDep;
Lon=reshape(Dataset(:,1),NumGrPt/NumDep,NumDep);
Lat=reshape(Dataset(:,2),NumGrPt/NumDep,NumDep);
Vs=reshape(Dataset(:,5),NumGrPt/NumDep,NumDep);


%% Interpolation
NumNode=length(GrLon)*length(GrLat);
[xq,yq]=meshgrid(GrLon,GrLat(end:-1:1));
VsMapTmp=zeros(NumNode,1);
VsMap=zeros(NumNode,length(Depth));
VsStr00=zeros(length(Dep),1);
for i=1:length(Dep)
    FVs=TriScatteredInterp(Lon(:,i),Lat(:,i),Vs(:,i),'natural');
    zq=FVs(xq,yq)';
    VsMapTmp(:,i)=reshape(zq,NumNode,1);
%     VsMapTmp(isnan(VsMapTmp(:,i))==1,i)=0;
    VsStr00(i,1)=mean(VsMapTmp(isnan(VsMapTmp(:,i))==0,i));
end
for i=1:NumNode
    VsMap(i,:)=interp1(Dep,VsMapTmp(i,:),Depth,'linear','extrap');
end
VsMean=interp1(Dep,VsStr00,Depth);
% VsMean=VsStr0;
OutCS(VsMap,VsMean,VphRslMdl(:,11),Grid,GrStep,0.01,Layer,strcat(PathOut,'PrevMdl/Cs/'));
OutVsMap(VsMean,VsMap,VsRslMdl,VphRslMdl,Layer,Grid,strcat(PathOut,'PrevMdl/Map/'));
