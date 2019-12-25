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


%% Body-wave seismic tomography
Fid=fopen('grid.result.txt','r');
Dataset=cell2mat(textscan(Fid,repmat('%f ',1,9),'HeaderLines',1));
Dep=unique(Dataset(:,3));
NumDepIn=length(Dep);
NumGrPt=length(Dataset(:,1));
NumNode1=NumGrPt/NumDepIn;
LonBst=reshape(Dataset(:,1),NumGrPt/NumDepIn,NumDepIn);
LatBst=reshape(Dataset(:,2),NumGrPt/NumDepIn,NumDepIn);
Vp=reshape(Dataset(:,4),NumGrPt/NumDepIn,NumDepIn);
Vs=reshape(Dataset(:,5),NumGrPt/NumDepIn,NumDepIn);
VpVs=reshape(Dataset(:,6),NumGrPt/NumDepIn,NumDepIn);
Res=reshape(Dataset(:,9),NumGrPt/NumDepIn,NumDepIn);

%% Surface-wave ambient noise tomography
load('RRT.mat');
Depth=[0;cumsum(Layer(1:end-1))];
NumNode2=size(Grid,1);
VsIntp=zeros(NumNode2,NumDepIn);
RslIntp=repmat(VphRslMdl(:,11),1,NumDepIn);
for i=1:NumNode2
    VsIntp(i,:)=interp1(Depth,VsStrMdl(i,:),Dep,'linear','extrap');
end

%% JMA Model
Fid=fopen('vjma2001','r');
Dataset=textscan(Fid,repmat('%f ',1,3));
fclose(Fid);
VpJMA=interp1(Dataset{3},Dataset{1},Dep);
VsJMA=interp1(Dataset{3},Dataset{2},Dep);

%% Aftershocks
% 3 month after the main shock (Determined by Okada et al., 2012)
Fid=fopen('hypo.txt','r');
Dataset=textscan(Fid,repmat('%f ',1,13),'HeaderLines',1);
fclose(Fid);
LatPst=Dataset{7};
LonPst=Dataset{8};
DepPst=Dataset{9};
MagPst=Dataset{10};
NumPost=length(LatPst);

%% Reflector
Fid=fopen('RefInfo.txt','r');
Info=textscan(Fid,repmat('%f ',1,9));
fclose(Fid);
RefData=cell2mat(Info);
NumRef=length(RefData);

%% Output
PathOutput=strcat(PathOut,'Cmp/');
if exist(PathOutput,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOutput,'s');
end
mkdir(PathOutput);

for i=1:NumDepIn
    Fid(1)=fopen(strcat(PathOutput,num2str(Dep(i),'%0.1f'),'km_Bst.txt'),'w');    
    Fid(2)=fopen(strcat(PathOutput,num2str(Dep(i),'%0.1f'),'km_Ant.txt'),'w');
    Fid(3)=fopen(strcat(PathOutput,num2str(Dep(i),'%0.1f'),'km_post.txt'),'w');
    Fid(4)=fopen(strcat(PathOutput,num2str(Dep(i),'%0.1f'),'km_ref.txt'),'w');
    VsPtb1=(Vs(:,i)/VsJMA(i)-1)*100;
    VsPtb2=(VsIntp(:,i)/VsJMA(i)-1)*100;
    for j=1:NumNode1
        fprintf(Fid(1),[repmat('%f ',1,4),'\n'],...
            LonBst(j,i),LatBst(j,i),Vs(j,i),VsPtb1(j));
    end
    for j=1:NumNode2
        fprintf(Fid(2),[repmat('%f ',1,5),'\n'],...
            Grid(j,:),VsIntp(j,i),VsPtb2(j),RslIntp(j,i));        
    end
    for j=1:NumPost
        if  LatPst(j) > min(Grid(:,2)) && LatPst(j) < max(Grid(:,2)) && ...
            LonPst(j) > min(Grid(:,1)) && LonPst(j) < max(Grid(:,1)) && ...
            abs(DepPst(j)-Dep(i)) < 0.5 && MagPst(j) > 0
            fprintf(Fid(3),[repmat('%.3f ',1,3),'%.1f\n'],...
                LonPst(j),LatPst(j),DepPst(j),MagPst(j));
        end
    end
    for j=1:NumRef
        if abs(RefData(j,5)-Dep(i)) < 0.5
            fprintf(Fid(4),[repmat('%.3f ',1,5),'\n'],RefData(j,[1,3,5,7,8]));
        end
    end
    fclose('all');
end




