%% Initial Setting
clear all;
close all;
addpath(genpath(pwd));
load('TMG.mat');


%% Produce File input for FMM tomography (by Rawlinson)
Fid(1)=fopen('source.dat','w');
Fid(2)=fopen('receivers.dat','w');
Fid(3)=fopen('otimes.dat','w');

FileId=fopen('Vertical.txt');
Channel=textscan(FileId,['%s%f%f%s%s%f%f%f%s',repmat('%f ',1,9)]);
fclose(FileId);
StLat=Channel{14};
StLon=Channel{15};
NumSt=length(StLat);
fprintf(Fid(1),'%d \n',NumSt);
fprintf(Fid(2),'%d \n',NumSt);
for i=1:NumSt
    fprintf(Fid(1),'%.2f %.2f \n',StLat(i),StLon(i));
    fprintf(Fid(2),'%.2f %.2f \n',StLat(i),StLon(i));
end

%% Identify station
RdErr=0.15;
Flag=zeros(NumSt,NumSt);
TrTime=zeros(NumSt,NumSt);
AprErr=zeros(NumSt,NumSt);
k=11; % Central frequency

NonZero=find(VphPthMdlIn(:,k)~=0, 1, 'last' );
for i=1:NonZero
    m=IdxStPair(i,1,k);
    n=IdxStPair(i,2,k);
     if m > n
        Flag(m,n)=1;
        TrTime(m,n)=StDist(i,k)/VphPthMdlIn(i,k);
        AprErr(m,n)=RdErr*StDist(i,k)/(VphPthMdlIn(i,k)^2-RdErr^2);        
    else
        Flag(n,m)=1;
        TrTime(n,m)=StDist(i,k)/VphPthMdlIn(i,k);
        AprErr(n,m)=RdErr*StDist(i,k)/(VphPthMdlIn(i,k)^2-RdErr^2);                
    end    
end

for m=1:NumSt
    for n=1:NumSt
        if Flag(m,n)==1
            fprintf(Fid(3),'%d %.4f %.2f \n',1,TrTime(m,n),AprErr(m,n));
        else
            fprintf(Fid(3),'%d %.4f %.2f \n',0,0,0.05);
        end        
    end
end


fclose('all');
