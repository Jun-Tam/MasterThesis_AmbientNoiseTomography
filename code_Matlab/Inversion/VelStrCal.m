function [VsCal,VsErr,VsRsl,VsFit,VsCmp,VphCal,VphErr,RmsItr,IdxOpt]=...
    VelStrCal(Vs0,VphObs,VphErr0,PriMdlErr,ItrStr,CrDist,Layer,Freq,PathK)

%% Initial Setting
C=[+0.9409;+2.0947;-0.8206;+0.2683;-0.0251;];
VpRho=1.741;
Robust=1.0;
IdxMinItr=10;
NumFreq=length(Freq);
NumLyr=length(Layer);
Depth=[0;cumsum(Layer(1:end-1))];
VsTmp=zeros(NumLyr,ItrStr);
VphTmp=zeros(NumFreq,ItrStr);
VsRslTmp=zeros(NumLyr,ItrStr);
VsErrTmp=zeros(NumLyr,ItrStr);
VphErrTmp=zeros(NumFreq,ItrStr);
RmsItr=zeros(1,ItrStr);
VsCmpTmp=zeros(1,ItrStr);

%% Inversion
for i=1:ItrStr
    %% Initial Model
    if i==1
        Vs0(end)=Vs0(end)*Robust;
        Vp=C(1)+C(2)*Vs0+C(3)*Vs0.^2+C(4)*Vs0.^3+C(5)*Vs0.^4;
        Rho=VpRho*Vp.^0.25;
        Fid=fopen(strcat('Model.k'),'w+');
        for j=1:NumLyr
            fprintf(Fid,[repmat('%6.2f ',1,4),'\n'],Layer(j),Rho(j),Vp(j),Vs0(j));
        end
        fclose(Fid);
        Vsi=Vs0;
    end
    
    %% Kernel Calculation
    eval(['!',PathK]);
    Fid(1)=fopen('RwPhVel.k');
    Fid(2)=fopen('RwPhSwKernel.k');
    Data1=textscan(Fid(1),'%f\n');
    Data2=textscan(Fid(2),[repmat('%f ',1,NumFreq),'\n']);
    fclose('all');
    VphMdl=Data1{1};
    K=zeros(NumFreq,NumLyr);
    for j=1:NumFreq
        KernelTmp=Data2{j};
        K(NumFreq-j+1,1:NumLyr)=KernelTmp(1+3*(0:(NumLyr-1)));
    end
    
    %% Inversion
    Ud=VphObs(end:-1:1);
    Um=VphMdl(end:-1:1);
    Cd=diag((VphErr0*1.0).^2);
    if CrDist==0 % Without smoothing
        Cm=PriMdlErr^2*eye(NumLyr);
    else         % With smoothing
        Cm=zeros(NumLyr);
        for l=1:NumLyr
            for m=1:NumLyr
                Denom=2*CrDist(l)*CrDist(m);
                Numer=-(Depth(l)-Depth(m))^2;
                Cm(l,m)=PriMdlErr^2*exp(Numer/Denom);
            end
        end
    end
    Vsi=Vs0+Cm*K'*pinv(Cd+K*Cm*K')*(Ud-Um+K*(Vsi-Vs0));
    CmPost=Cm-Cm*K'*pinv(Cd+K*Cm*K')*K*Cm;
    Rsl=eye(NumLyr)-CmPost*pinv(Cm);
    VsTmp(:,i)=Vsi;
    VphTmp(:,i)=Um(end:-1:1);
    VsRslTmp(:,i)=diag(Rsl);
    VsErrTmp(:,i)=sqrt(diag(CmPost));
    VphErrTmp(:,i)=sqrt(diag(K*Cm*K'));
    RmsItr(i)=rms(Ud-Um);
    VsCmpTmp(i)=var(Vsi-Vs0);
    
    %% Output
    
    fid=fopen(strcat('Model.k'),'w+');
    Vsi(end)=Vsi(end)*Robust;    
    Vp=C(1)+C(2)*Vsi+C(3)*Vsi.^2+C(4)*Vsi.^3+C(5)*Vsi.^4;
    Rho=VpRho*Vp.^0.25;
    for j=1:NumLyr
        fprintf(fid,[repmat('%6.2f ',1,4),'\n'],Layer(j),Rho(j),Vp(j),Vsi(j));
    end
    fclose(fid);    
end

OutKernel(Depth,Freq,K);
[~,IdxMin]=min(RmsItr(IdxMinItr:end));
IdxOpt=IdxMin+IdxMinItr-1;
% IdxOpt=ItrStr;
VsCal=VsTmp(:,IdxOpt);
VphCal=VphTmp(:,IdxOpt);
VsRsl=VsRslTmp(:,IdxOpt);
VsErr=VsErrTmp(:,IdxOpt);
VphErr=VphErrTmp(:,IdxOpt);
VsFit=RmsItr(IdxOpt);
VsCmp=VsCmpTmp(IdxOpt);
delete('*.k');
delete('fort*');

end

