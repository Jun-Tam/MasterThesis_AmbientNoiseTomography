function Resolution=TmgInvResl(R,Grid,GrStep)
% Resolution : Resolution map

NumNode=sqrt(length(R));
R=reshape(R,NumNode,NumNode);
RMapTmp=zeros(1,NumNode);
NumRad=10;
NumAzm=21;
RadMax=[.5:.01:1, 1.1:.1:10];
NumRadGrd=length(RadMax);

for i=1:NumNode
    if nnz(R(:,i))~=0
        Misfit=zeros(NumRadGrd,1);
        [Hobs,GrLon,GrLat]=Grid2Mesh(R(:,i),Grid,GrStep);
        [Height,IdxSummit]=max(R(:,i),[],1);
        PosSummit=Grid(IdxSummit,:);
        Hobs=1-Hobs/max(Height);
        %% Grid search
        for j=1:NumRadGrd
            [Rad,Azm]=meshgrid(linspace(0,RadMax(j),NumRad),linspace(0,2*pi,NumAzm));
            Lon=PosSummit(1)+Rad.*cos(Azm);     % Position for data
            Lat=PosSummit(2)+Rad.*sin(Azm);     % Position for data
            Hmdl=1-Rad./RadMax(j);
            [LonMesh,LatMesh]=meshgrid(GrLon,GrLat);
            Hobs2=reshape(interp2(LonMesh,LatMesh,Hobs,Lon,Lat),numel(Lon),1);
            Hmdl2=reshape(Hmdl,numel(Lon),1);
            IdxNaN=isnan(Hobs2)==1;
            Misfit(j)=rms(Hmdl2(IdxNaN~=1)-Hobs2(IdxNaN~=1));
        end
        Misfit=moving(Misfit,3);
        RMapTmp(i)=RadMax(Misfit==min(Misfit));
    end
end
LatDis=distance(min(GrLat),mean(GrLon),max(GrLat),mean(GrLon));
LonDis=distance(mean(GrLat),min(GrLon),mean(GrLat),max(GrLon));
CoefRad=(LonDis/length(GrLon)+LatDis/length(GrLat))/2/1000/360*(2*pi)*earthRadius;
Resolution=CoefRad*reshape(RMapTmp,NumNode,1);
IdxRmv=Resolution==0;
Resolution(IdxRmv)=NaN;

%% Check
% figure;
% [Image,~,~]=Grid2Mesh(Resolution,Grid,GrStep);
% imagesc(GrLon,GrLat,Image);
% set(gca,'ydir','normal','clim',[5 10]);
%         figure;
%         subplot(1,2,1);
%         surf(GrLon,GrLat,Hobs);
%         set(gca,'zlim',[0 1]);
%         subplot(1,2,2);
%         surf(GrLon,GrLat,Hmdl);
%         set(gca,'zlim',[0 1]);

% warning('off','MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId')
% [Xmesh,Ymesh]=meshgrid(GrLon,GrLat);
% X=reshape(Lon,numel(Lon),1);
% Y=reshape(Lat,numel(Lat),1);
% Z=reshape(Rad/RMapTmp(i),numel(Rad),1);
% F=TriScatteredInterp(X,Y,Z,'natural');
% Hmdl=F(Xmesh,Ymesh);
% Misfit(i)=rms(rms(Hmdl-Hobs));
end

