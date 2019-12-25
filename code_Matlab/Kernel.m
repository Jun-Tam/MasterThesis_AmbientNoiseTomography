%% Initial Setting
Layer=ones(11,1);
Depth=0:1:10;
Freq=0.15:0.025:0.40;
Vs=[2.844;3.012;3.157;3.278;3.375;3.431;3.451;3.471;3.491;3.511;3.531;];
C=[+0.9409;+2.0947;-0.8206;+0.2683;-0.0251;];
VpRho=1.741;
Vp=C(1)+C(2)*Vs+C(3)*Vs.^2+C(4)*Vs.^3+C(5)*Vs.^4;
Rho=VpRho*Vp.^0.25;
NumLyr=length(Layer);
NumFrq=length(Freq);
PathK=strcat(pwd,'/RayleighM');
PathOutput=strcat(pwd,'/ANT/Output/Kernel/');
setenv('DYLD_LIBRARY_PATH','/usr/local/bin/');

%% Calculation
Fid=fopen(strcat('Model.k'),'w+');
for j=1:NumLyr
    fprintf(Fid,[repmat('%6.2f ',1,4),'\n'],Layer(j),Rho(j),Vp(j),Vs(j));
end
fclose(Fid);
eval(['!',PathK]);
Fid(1)=fopen('RwPhVel.k');
Fid(2)=fopen('RwPhSwKernel.k');
Data1=textscan(Fid(1),'%f\n');
Data2=textscan(Fid(2),[repmat('%f ',1,NumFrq),'\n']);
fclose('all');
VphMdl=Data1{1};
K=zeros(NumFrq,NumLyr);
for j=1:NumFrq
    KernelTmp=Data2{j};
    K(NumFrq-j+1,1:NumLyr)=KernelTmp(1+3*(0:(NumLyr-1)));
end

%% Output
for i=1:NumFrq
    Fid=fopen(strcat(PathOutput,num2str(Freq(i),'%2.3f'),'Hz.txt'),'w');
    for j=1:NumLyr
        fprintf(Fid,'%.3f %.3f\n',K(NumFrq-i+1,j),Depth(j));
    end
    fclose(Fid);
end
