function OutDisp(VphMsr,VphRef,VrRef,VphPthMdlIn,VphMsrVr,Pos,Freq,PathOut)

NumFrq=length(Freq);
NumData=size(VphMsr,2);
if exist(PathOut,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOut,'s');
end
mkdir(PathOut);

%% Dispersion curves (Complete path)
for i=1:NumData
    if sum(VphMsr(:,i)~=0)==NumFrq
        St1Pos=strcat(num2str(Pos(i,1),'%.4f'),'_',num2str(Pos(i,2),'%.4f'));
        St2Pos=strcat(num2str(Pos(i,3),'%.4f'),'_',num2str(Pos(i,4),'%.4f'));        
        FileName=strcat(PathOut,St1Pos,'_',St2Pos,'.txt');
        Fid=fopen(FileName,'w');
        for j=1:NumFrq
            fprintf(Fid,'%.3f %.3f %.3f\n',Freq(j),VphMsr(j,i),VphMsrVr(j,i));
        end
        fclose(Fid);
    end
end

%% JMA Model
Layer=[ones(10,1);2;];
NumLyr=length(Layer);
Depth=[0;cumsum(Layer(1:NumLyr-1))];
Fid(1)=fopen('vjma2001','r');
PathR=strcat(pwd,'/RayleighJMA');
Dataset=textscan(Fid(1),repmat('%f ',1,3));
VpJMA=Dataset{1};
VsJMA=Dataset{2};
DepJMA=Dataset{3};
Vp=interp1(DepJMA,VpJMA,Depth);
Vs=interp1(DepJMA,VsJMA,Depth);
Rho=1.741*Vp.^0.25;
FreqJMA=0.15:0.025:0.50;
NumFrqJMA=length(FreqJMA);
Fid(2)=fopen(strcat('Model.k'),'w+');
for i=1:NumLyr
    fprintf(Fid(2),[repmat('%6.2f ',1,4),'\n'],Layer(i),Rho(i),Vp(i),Vs(i));
end
eval(['!',PathR]);
Fid(3)=fopen('RwPhVel.k');
Fid(4)=fopen(strcat(PathOut,'DispJMA.txt'),'w');
Data1=textscan(Fid(3),'%f\n');
VphJMA=Data1{1};

% FidGr=fopen('RwGrVel.k');
% DataGr=textscan(FidGr,'%f\n');
% VphJmaGr=DataGr{1};
for i=1:NumFrqJMA
    fprintf(Fid(4),'%.3f %.3f\n',FreqJMA(i),VphJMA(i));
end
delete('./*.k');

%% S-wave velocity (JMA2001)
Fid(5)=fopen(strcat(PathOut,'VsJMA.txt'),'w');
for i=1:NumLyr
    fprintf(Fid(5),'%.3f %.3f %.3f %.3f\n',Depth(i)-0.5,Vp(i),Vs(i),Rho(i));
    fprintf(Fid(5),'%.3f %.3f %.3f %.3f\n',Depth(i)+0.5,Vp(i),Vs(i),Rho(i));    
end
% for i=1:NumLyr
%     fprintf(Fid,'%.3f %.3f %.3f %.3f\n',Depth(i),Vp(i),Vs(i),Rho(i));
% end

%% Dispersion curve (Reference, Average, JMA2001)
Fid(6)=fopen(strcat(PathOut,'DispRef.txt'),'w');
Fid(7)=fopen(strcat(PathOut,'DispAvg.txt'),'w');
Fid(8)=fopen(strcat(PathOut,'DispJmaShort.txt'),'w');
for i=1:NumFrq
    NumPath=find(VphPthMdlIn(:,i)~=0, 1, 'last' );
    VphAvg=mean(VphPthMdlIn(1:NumPath,i));
    VphStd=std(VphPthMdlIn(1:NumPath,i));
    fprintf(Fid(6),'%.3f %.3f %.3f\n',Freq(i),VphRef(i),VrRef(i));    
    fprintf(Fid(7),'%.3f %.3f %.3f\n',Freq(i),VphAvg,VphStd);
    fprintf(Fid(8),'%.3f %.3f\n',Freq(i),VphJMA(i));    
end
fclose('all');

end
