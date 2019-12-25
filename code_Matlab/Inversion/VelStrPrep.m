function VsModel=VelStrPrep(Layer,Freq,VphRef,VphErrIn,CrDist,DmpStr,...
    ItrStr,FlagVs0,PathK)

Depth=[0;cumsum(Layer(1:end-1))];

%% Initial Model Calculation
if     FlagVs0==1 % Model 1 (Reference Dispersion curve)
    Fid0=fopen('vjma2001','r');
    Dataset=textscan(Fid0,repmat('%f ',1,3));
    fclose(Fid0);
    VsJMA=Dataset{2};
    DepJMA=Dataset{3};
    VsModel0=interp1(DepJMA,VsJMA,Depth);
    [VsModel,~,~,~,~,~,~,~,~]=VelStrCtrl(VsModel0,VphRef',...
        VphErrIn,DmpStr,ItrStr,CrDist,Layer,Freq,PathK);
elseif FlagVs0==2 % Model 2 (JMA)
    Fid0=fopen('vjma2001','r');
    Dataset=textscan(Fid0,repmat('%f ',1,3));
    fclose(Fid0);
    VsJMA=Dataset{2};
    DepJMA=Dataset{3};
    VsModel=interp1(DepJMA,VsJMA,Depth)';
end

end