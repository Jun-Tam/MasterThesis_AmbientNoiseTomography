function [VphMapRcvIn,VsStrRcvIn]=RestoreRcv(VsStr,Layer,Freq,Grid,PathK)

%% Setting
C=[+0.9409;+2.0947;-0.8206;+0.2683;-0.0251;];
VpRho=1.741;
NumNode=size(Grid,1);
NumLyr=length(Layer);
NumFreq=length(Freq);
VphMapRcvIn=zeros(NumNode,NumFreq);
VsStrRcvIn=VsStr;

%% Synthetic model
for i=1:NumNode
    Fid=fopen(strcat('Model.k'),'w+');
    Vs=VsStr(i,:);
    Vp=C(1)+C(2)*Vs+C(3)*Vs.^2+C(4)*Vs.^3+C(5)*Vs.^4;
    Rho=VpRho*Vp.^0.25;
    for j=1:NumLyr
        fprintf(Fid,[repmat('%6.2f ',1,4),'\n'],Layer(j),Rho(j),Vp(j),Vs(j));
    end
    eval(['!',PathK]);
    Fid=fopen('RwPhVel.k');
    Data1=textscan(Fid,'%f\n');
    fclose('all');
    VphMapRcvIn(i,:)=Data1{1};
    delete('./*.k');
end

end

