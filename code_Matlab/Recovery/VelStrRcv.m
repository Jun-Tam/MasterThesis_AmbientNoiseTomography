function [VphPthRcvIn,VphMapRcvIn,VsStrRcvIn]=VelStrRcv(IdxDep,StDist,...
    StAzim,StPos,VsStr,Grid,GrSize,Lon,Lat,Layer,Freq,NumSect,NumBinAzm,...
    FlagSect,PathK,FlagRcv,Noise)

tic;
if FlagRcv==1      % Spike Resolution Test
    [VphMapRcvIn,VsStrRcvIn]=SpikeRcv(VsStr,Layer,Freq,Grid,IdxDep,PathK);
elseif FlagRcv==2  % Restoring Resolution Test
    [VphMapRcvIn,VsStrRcvIn]=RestoreRcv(VsStr,Layer,Freq,Grid,PathK);
end

%% Path processing
NumFreq=length(Freq);
VphPthRcvIn=zeros(1000,NumFreq);
if numel(GrSize)==1     % Regular grid
    [GrStep,Grid,IdxCorner,NumWgtPos]=RglGrid(Lon,Lat,GrSize);
else                    % Adaptive grid
    [GrStep,Grid,IdxCorner,NumWgtPos]=AdpGrid(Lon,Lat,GrSize);    
end
for i=1:NumFreq
    IdxNonZero=StDist(:,i)~=0;
    NumData=sum(IdxNonZero);
    Dist0=StDist(IdxNonZero,i);
    Azim0=StAzim(IdxNonZero,i);
    Pos0=StPos(IdxNonZero,1:4,i);    
    [Weight,IdxPosWgt,IdxPosBox,NumDiv,~,~,~]=TmgInvPath(Pos0,Dist0,...
        Azim0,IdxCorner,GrStep,Grid,NumWgtPos,NumSect,NumBinAzm,FlagSect);
    for j=1:NumData
        tmp=0;
        for k=1:NumDiv(j)
            kmax=NumWgtPos(IdxPosBox(k,j));
            dPath=Dist0(j)/NumDiv(j);
            dVel=sum(VphMapRcvIn(IdxPosWgt(k,j,1:kmax),i).*...
                reshape(Weight(k,j,1:kmax),kmax,1));
            tmp=tmp+dPath/dVel;
        end
        VphPthRcvIn(j,i)=Dist0(j)/(tmp+Noise*randn(1));
    end
end
toc;

end
