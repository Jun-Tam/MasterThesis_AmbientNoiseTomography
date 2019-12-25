function [VphMapRcvIn,VsStrRcvIn] = SpikeRcv(VsStr,Layer,Freq,Grid,IdxDep,PathK)

%% Setting
if strcmp(IdxDep,'24');
    LyrLow=2:4;
elseif strcmp(IdxDep,'46')==1;
    LyrLow=4:6;
elseif strcmp(IdxDep,'68')==1;
    LyrLow=6:8;
end
LonLow=[140.67,140.73];
LatLow=[38.77,38.83];
AmpLow=20;
C=[+0.9409;+2.0947;-0.8206;+0.2683;-0.0251;];
VpRho=1.741;
NumNode=size(Grid,1);
NumLyr=length(Layer);
NumFreq=length(Freq);
Coeff=ones(1,numel(Layer));
Coeff(LyrLow)=1-AmpLow/100;
VsRcv=zeros(2,NumLyr);
VphRcv=zeros(2,NumFreq);
VphMapRcvIn=zeros(NumNode,NumFreq);
VsStrRcvIn=zeros(NumNode,NumLyr);
VsRcv(1,:)=mean(VsStr,1);
VsRcv(2,:)=mean(VsStr,1).*Coeff;

%% Synthetic model
for i=1:2
    Fid=fopen(strcat('Model.k'),'w+');
    Vs=VsRcv(i,:);
    Vp=C(1)+C(2)*Vs+C(3)*Vs.^2+C(4)*Vs.^3+C(5)*Vs.^4;
    Rho=VpRho*Vp.^0.25;
    for j=1:NumLyr
        fprintf(Fid,[repmat('%6.2f ',1,4),'\n'],Layer(j),Rho(j),Vp(j),Vs(j));
    end
    eval(['!',PathK]);
    Fid=fopen('RwPhVel.k');
    Data1=textscan(Fid,'%f\n');
    fclose('all');
    VphRcv(i,:)=Data1{1};
    delete('./*.k');
end
for i=1:NumNode
    if  Grid(i,1) >= LonLow(1) && Grid(i,1) <= LonLow(2) && ...
        Grid(i,2) >= LatLow(1) && Grid(i,2) <= LatLow(2)
        VphMapRcvIn(i,:)=VphRcv(2,:);
        VsStrRcvIn(i,:)=VsRcv(2,:);
    else
        VphMapRcvIn(i,:)=VphRcv(1,:);
        VsStrRcvIn(i,:)=VsRcv(1,:);
    end
end

end

