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

%% Input
Fid(1)=fopen('Vs3100.txt','r');
Fid(2)=fopen('MeshCode','r');
DepVs3100=cell2mat(textscan(Fid(1),repmat('%d ',1)));
MeshCode=textscan(Fid(2),repmat('%s ',1));
fclose('all');
Code=cellstr(char(MeshCode{1}));
NumNode=length(Code);
Lon=zeros(NumNode,1);
Lat=zeros(NumNode,1);
for i=1:NumNode
    Str=char(Code(i));
    Lat1=str2double(Str(1:2))/1.5;  % Longitude
    Lon1=str2double(Str(3:4))+100;  % Latitude
    Lat2=str2double(Str(5))*5;
    Lon2=str2double(Str(6))*7.5;
    Lat3=str2double(Str(7))*30;
    Lon3=str2double(Str(8))*45;
    Lat(i)=Lat1+Lat2/60+Lat3/3600;
    Lon(i)=Lon1+Lon2/60+Lon3/3600;
end

%% Output
Fid=fopen('BedRockDepth.txt','w');
for i=1:NumNode
    fprintf(Fid,'%.4f %.4f %d \n',Lon(i),Lat(i),DepVs3100(i));
end
fclose(Fid);
