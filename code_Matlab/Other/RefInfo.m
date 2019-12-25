%% Initial Setting
clear all;
close all;
Path='/Users/Jun/Dropbox';
PathHome=strcat(Path,'/ANT');
PathIn=strcat(PathHome,'/Input/');
PathOut=strcat(PathHome,'/Output/');
addpath(genpath(pwd));
addpath(genpath(strcat(Path,'/Code/Matlab')));
addpath(genpath(strcat(Path,'/InfoSource')));
setenv('DYLD_LIBRARY_PATH','/usr/local/bin/');
cd(PathHome);

Layer=[ones(8,1);2;];
Depth=[0;cumsum(Layer(1:end-1))];
Lat=[38.54, 39.10];
Lon=[140.36,141.08];

%% Input
i=1;
Fid=fopen('/Users/Jun/Dropbox/ANT/Reference/Reflector.txt','r');
while fgets(Fid)~=-1
    Info=textscan(Fid,repmat('%f ',1,5));
    if isempty(cell2mat(Info)) == 0
        InfoMat=cell2mat(Info);        
        NumRow=size(InfoMat,1)-1;
        DataRef(i:i+NumRow-1,:)=InfoMat(1:NumRow,:);
        i=i+NumRow;
    else
        fgets(Fid);
    end
end
fclose(Fid);

%% Selection
IdxUse=DataRef(:,1) > min(Lon)   & DataRef(:,1) < max(Lon) & ...
       DataRef(:,2) > min(Lat)   & DataRef(:,2) < max(Lat) & ...
       DataRef(:,3) > min(Depth) & DataRef(:,3) < max(Depth);
RefLon=DataRef(IdxUse,1);
RefLat=DataRef(IdxUse,2);
RefDep=DataRef(IdxUse,3);
RefStr=DataRef(IdxUse,4);
RefDip=DataRef(IdxUse,5);

%% Statistical processing
RefStrTmp=RefStr(1);
RefDipTmp=RefDip(1);
RefLonTmp=zeros(100,1);
RefLatTmp=zeros(100,1);
RefDepTmp=zeros(100,1);
RefLonMean=zeros(100,1);
RefLatMean=zeros(100,1);
RefDepMean=zeros(100,1);
RefLonVar=zeros(100,1);
RefLatVar=zeros(100,1);
RefDepVar=zeros(100,1);
RefStrRef=zeros(100,1);
RefDipRef=zeros(100,1);
RefCount=zeros(100,1);
j=0;
k=0;
for i=2:sum(IdxUse)
    if RefStr(i)==RefStrTmp && RefDip(i)==RefDipTmp
        k=k+1;
        RefLonTmp(k)=RefLon(i);
        RefLatTmp(k)=RefLat(i);
        RefDepTmp(k)=RefDep(i);
    else
        j=j+1;
        IdxZero=RefLonTmp==0;
        RefLonTmp(IdxZero)=[];
        RefLatTmp(IdxZero)=[];
        RefDepTmp(IdxZero)=[];
        RefLonMean(j)=mean(RefLonTmp);
        RefLatMean(j)=mean(RefLatTmp);
        RefDepMean(j)=mean(RefDepTmp);
        RefLonVar(j)=var(RefLonTmp);
        RefLatVar(j)=var(RefLatTmp);
        RefDepVar(j)=var(RefDepTmp);
        RefStrRef(j)=RefStrTmp;
        RefDipRef(j)=RefDipTmp;
        RefCount(j)=k;
        RefStrTmp=RefStr(i);
        RefDipTmp=RefDip(i);
        k=0;
        RefLonTmp=zeros(100,1);
        RefLatTmp=zeros(100,1);
        RefDepTmp=zeros(100,1);
    end
end
NumRef=j;
IdxZero=RefLonMean==0;
RefLonMean(IdxZero)=[];
RefLatMean(IdxZero)=[];
RefDepMean(IdxZero)=[];
RefLonVar(IdxZero)=[];
RefLatVar(IdxZero)=[];
RefDepVar(IdxZero)=[];
RefStrRef(IdxZero)=[];
RefDipRef(IdxZero)=[];
RefCount(IdxZero)=[];

%% Output
Fid=fopen(strcat('RefInfo.txt'),'w');
for i=1:NumRef
       fprintf(Fid,[repmat('%f ',1,8),'%d ','\n'],...
           RefLonMean(i),RefLonVar(i),RefLatMean(i),RefLatVar(i),...
           RefDepMean(i),RefDepVar(i),RefStrRef(i),RefDipRef(i),RefCount(i));    
end
fclose('all');
