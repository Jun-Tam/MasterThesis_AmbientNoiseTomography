function OutSPAC(Dist,Azim,VR,VphVR,FreqVR,AmpVR,SpReal,SpImag,...
    Freq,VphMsr,VelGrRange,VelGrStep,Cmp,PathOutput)

if exist(PathOutput,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOutput,'s');
end
mkdir(PathOutput);

DistPlot=0:0.1:60;
NumPlot=length(DistPlot);
NumData=length(Dist);
NumFreqVR=length(FreqVR);
VR(VR < 1e-3)=0.001;


%% Fitting
[x,y]=meshgrid(FreqVR,VelGrRange(1):VelGrStep:VelGrRange(2));
NumGr=numel(VR);
x=reshape(x',NumGr,1);
y=reshape(y',NumGr,1);
z=reshape(VR,NumGr,1);
Fid(1)=fopen(strcat(PathOutput,'/VarR.txt'),'w');
Fid(2)=fopen(strcat(PathOutput,'/Ref.txt'),'w');
for i=1:NumGr
    fprintf(Fid(1),[repmat('%.3f ',1,3),'\n'],x(i),y(i),z(i));
end
for i=1:NumFreqVR
    fprintf(Fid(2),[repmat('%.3f ',1,2),'\n'],FreqVR(i),VphVR(i));
end
fclose('all');

%% Dispersion curve
for i=1:size(VphMsr,1)
    Fid=fopen(strcat(PathOutput,'/DspCrv',num2str(i,'%d'),'.txt'),'w');
    for j=1:length(Freq)
        fprintf(Fid,[repmat('%0.3f ',1,2),'\n']',Freq(j),VphMsr(i,j));
    end
    fclose('all');
end

%% Spatial spectrum distribution
for i=1:NumFreqVR
    Nrm=max(abs(SpReal(:,i)));
    SpRealNrm=SpReal(:,i)/Nrm;
    SpImagNrm=SpImag(:,i)/Nrm;
    WnDist=2*pi*FreqVR(i)/VphVR(i)*DistPlot;
    if strcmp(Cmp,'ZZ')==1
        ApprxPlot=real(AmpVR(i))*besselj(0,WnDist)/Nrm;
    elseif strcmp(Cmp,'TT')==1 || strcmp(Cmp,'RR')==1
        ApprxPlot=AmpVR(i)*((besselj(0,WnDist)-besselj(2,WnDist))/2)/Nrm;
    end
    Path=strcat(PathOutput,'/',num2str(FreqVR(i),'%2.3f'));
    Fid(1)=fopen(strcat(Path,'Hz(Real)Sp.txt'),'w');
    Fid(2)=fopen(strcat(Path,'Hz(Imag)Sp.txt'),'w');
    Fid(3)=fopen(strcat(Path,'Hz(Real)Ds.txt'),'w');
    Fid(4)=fopen(strcat(Path,'Hz(Imag)Ds.txt'),'w');    
    Fid(5)=fopen(strcat(Path,'Hz(Real)Apprx.txt'),'w');
    for j=1:NumData
        fprintf(Fid(1),[repmat('%.3f ',1,3),'\n'],...
            Dist(j).*cosd(Azim(j)),Dist(j).*sind(Azim(j)),SpRealNrm(j));
        fprintf(Fid(1),[repmat('%.3f ',1,3),'\n'],...
            -Dist(j).*cosd(Azim(j)),-Dist(j).*sind(Azim(j)),SpRealNrm(j+NumData));
        fprintf(Fid(2),[repmat('%.3f ',1,3),'\n'],...
            Dist(j).*cosd(Azim(j)),Dist(j).*sind(Azim(j)),SpImagNrm(j));
        fprintf(Fid(2),[repmat('%.3f ',1,3),'\n'],...
            -Dist(j).*cosd(Azim(j)),-Dist(j).*sind(Azim(j)),SpImagNrm(j+NumData));        
        fprintf(Fid(3),[repmat('%.3f ',1,2),'\n'],Dist(j),SpRealNrm(j));
        fprintf(Fid(3),[repmat('%.3f ',1,2),'\n'],Dist(j),SpRealNrm(j+NumData));
        fprintf(Fid(4),[repmat('%.3f ',1,2),'\n'],Dist(j),SpImagNrm(j));
        fprintf(Fid(4),[repmat('%.3f ',1,2),'\n'],Dist(j),SpImagNrm(j+NumData));        
    end
    for j=1:NumPlot
       fprintf(Fid(5),'%0.3f %0.3f\n',DistPlot(j),ApprxPlot(j)) ;
    end    
    fclose('all');
end

end
