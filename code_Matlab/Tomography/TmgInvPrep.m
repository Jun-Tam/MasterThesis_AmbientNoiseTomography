function [VphPthIn,VphMapOut,VphPthOut,VphMapIn]=TmgInvPrep(VphPthIn,VphRef,...
    Dist,Grid,GrStep,NumWgtPos,Weight,IdxPosWgt,IdxPosBox,NumDiv,FlagMdl)

%% Initialization
NumNode=size(Grid,1);
NumPath=length(VphPthIn);
VphPthOut=zeros(NumPath,1);

%% Initial Model
if FlagMdl==0; %VphPthIn;
    VphMapIn=zeros(NumNode,1);
else
    LowVelRate=30;
    Noise=0.2; %% Noise 0.2 s
    if FlagMdl==1 || FlagMdl==2 % CBT 0.04 by 0.04 or 0.08 by 0.08 
        NumAnm=FlagMdl;
        [~,GrLon,GrLat]=Grid2Mesh(Grid(:,1),Grid,GrStep);
        LngLon=length(GrLon);
        LngLat=length(GrLat);
        RndLon=round(LngLon/2)/NumAnm;
        RndLat=round(LngLat/2)/NumAnm;
        ChkBrd=checkerboard(NumAnm,RndLon,RndLat) > 0.5;
        ChkBrd=double(ChkBrd);
        if mod(LngLon,2)==1
            ChkBrd(:,1)=[];
        end
        if mod(LngLat,2)==1
            ChkBrd(1,:)=[];
        end
        VphMapIn=zeros(LngLon,LngLat);
        for n=1:LngLon
            for m=1:LngLat
                if ChkBrd(n,m)==1
                    VphMapIn(n,m)=VphRef*(1+LowVelRate/100);
                else
                    VphMapIn(n,m)=VphRef*(1-LowVelRate/100);
                end
            end
        end
        VphMapIn=Mesh2Grid(VphMapIn',Grid,GrStep);
    else % Spike Test
        RangeLon=[140.72,140.74]; % Naruko
        RangeLat=[38.72,38.74];
        Idx1=Grid(:,1) >= min(RangeLon) & Grid(:,1) <= max(RangeLon) & ...
             Grid(:,2) >= min(RangeLat) & Grid(:,2) <= max(RangeLat);
        RangeLon=[140.66,140.72]; % Onikobe1
        RangeLat=[38.80,38.83];
        Idx2=Grid(:,1) >= min(RangeLon) & Grid(:,1) <= max(RangeLon) & ...
             Grid(:,2) >= min(RangeLat) & Grid(:,2) <= max(RangeLat);
        Idx=Idx1 | Idx2;
        VphMapIn(1:NumNode,1)=VphRef;
        VphMapIn(Idx,1)=VphRef*(1-LowVelRate/100);
    end
    for i=1:NumPath
        tmp=0;
        for j=1:NumDiv(i)
            kmax=NumWgtPos(IdxPosBox(j,i));
            dPath=Dist(i)/NumDiv(i);
            dVel=sum(VphMapIn(IdxPosWgt(j,i,1:kmax)).*Weight(j,i,1:kmax));
            tmp=tmp+dPath/dVel;
        end
        VphPthIn(i)=Dist(i)/(tmp+Noise*randn(1));
    end
end
VphMapOut=repmat(VphRef,1,NumNode); % Homogeneous Initial Model

%% Model velocity on used path
for i=1:NumPath
    tmp=0;
    for j=1:NumDiv(i)
        kmax=NumWgtPos(IdxPosBox(j,i));
        dPath=Dist(i)/NumDiv(i);
        tmp=tmp+dPath/sum(VphMapOut(IdxPosWgt(j,i,1:kmax)).*Weight(j,i,1:kmax));
    end
    VphPthOut(i,1)=Dist(i)/tmp;
end

end

