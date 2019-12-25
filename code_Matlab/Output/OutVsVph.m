function OutVsVph(Grid,GrStep,Layer,Freq,VsStr,VsErr,VphMap,VphStr,VphErr,...
    VsStr0,VphRef,IdxDmp,Cmp,Fit,Dmp,FlagDmp,ItrStr,IdxItr,RmsItr,PathOutput)

%% Altitude Input
Fid=fopen('alt.txt','r');
Data=textscan(Fid,'%f %f %f\n');
AltLon=Data{1};
AltLat=Data{2};
Altitude=Data{3};
GrLonDsp=min(Grid(:,1)):min(GrStep):max(Grid(:,1));
GrLatDsp=min(Grid(:,2)):min(GrStep):max(Grid(:,2));
[xq,yq]=meshgrid(GrLonDsp,GrLatDsp);
F=TriScatteredInterp(AltLon,AltLat,Altitude,'natural'); 
AltMap=F(xq,yq);
AltNode=Mesh2Grid(AltMap,Grid,GrStep);

%% ------------------------------------------------------------------------
if exist(PathOutput,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOutput,'s');
end
mkdir(PathOutput);

NumDmp=length(Dmp);
NumNode=size(Grid,1);
NumFrq=length(Freq);
NumLyr=length(Layer);
Depth0=[0;cumsum(Layer(1:end-1))];

%% Vs profile
for i=1:NumNode
    FileName=strcat(num2str(Grid(i,1),'%0.2f'),'_',num2str(Grid(i,2),'%0.2f'));
    Fid(1)=fopen(strcat(PathOutput,FileName,'_Vs.txt'),'w');
    Fid(2)=fopen(strcat(PathOutput,FileName,'_VsErr.txt'),'w');    
    Fid(3)=fopen(strcat(PathOutput,FileName,'_Vph.txt'),'w');
    Fid(4)=fopen(strcat(PathOutput,FileName,'_Dmp.txt'),'w');
    Fid(5)=fopen(strcat(PathOutput,FileName,'_Itr.txt'),'w');
%     if isnan(AltNode(i))==0
%         Depth=Depth0-AltNode(i)/1000;
%     else
%         Depth=Depth0;
%     end
%     for j=1:NumLyr
%         fprintf(Fid(1),[repmat('%0.3f ',1,2),'\n'],Depth0(j)-0.5,VsStr(i,j));
%         fprintf(Fid(1),[repmat('%0.3f ',1,2),'\n'],Depth0(j)+0.5,VsStr(i,j));
%         fprintf(Fid(2),[repmat('%0.3f ',1,3),'\n'],Depth0(j),VsStr(i,j),VsErr(i,j));
%     end
    %% Vs profile
    for j=1:NumLyr
        fprintf(Fid(1),[repmat('%0.3f ',1,2),'\n'],Depth0(j),VsStr(i,j));
        fprintf(Fid(2),[repmat('%0.3f ',1,3),'\n'],Depth0(j),VsStr(i,j),VsErr(i,j));
    end
    
    %% Dispersion curve
    for j=1:NumFrq
        fprintf(Fid(3),[repmat('%0.3f ',1,4),'\n'],...
            Freq(j),VphMap(i,j),VphStr(i,j),VphErr(i,j));        
    end
        
    %% A priori model error (trade-off)
    k=IdxDmp(i);
    fprintf(Fid(4),'%.3f %f %f \n',Dmp(k),Fit(i,k),Cmp(i,k));
    for j=1:NumDmp
        if Cmp(i,j)~=0 && FlagDmp(i,j)==1
            fprintf(Fid(4),'%f %f %f \n',Dmp(j),Fit(i,j),Cmp(i,j));
        end
    end
    
    %% Iteration and Rms residual
    fprintf(Fid(5),'%d %f \n',IdxItr(i),RmsItr(i,IdxItr(i)));
    for j=1:ItrStr
        fprintf(Fid(5),'%d %f \n',j,RmsItr(i,j));
    end    
    fclose('all');
end

%% Reference for Vs profile and Dispersion curve
Fid(6)=fopen(strcat(PathOutput,'VsRef.txt'),'w');
Fid(7)=fopen(strcat(PathOutput,'VphRef.txt'),'w');
% for j=1:NumLyr-1
%     fprintf(Fid(5),[repmat('%0.3f ',1,2),'\n'],Depth0(j),VsStr0(j));
%     fprintf(Fid(5),[repmat('%0.3f ',1,2),'\n'],Depth0(j+1),VsStr0(j));
% end
for j=1:NumLyr
    fprintf(Fid(6),[repmat('%0.3f ',1,2),'\n'],Depth0(j),VsStr0(j));
end
for j=1:NumFrq
    fprintf(Fid(7),[repmat('%0.3f ',1,2),'\n'],Freq(j),VphRef(j));
end
fclose('all');

end
