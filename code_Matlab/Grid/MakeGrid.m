function Grid = MakeGrid(Lon,Lat,GrSize)

%% Initial Setting
Eps=1e-4;
Grid=zeros(10000,2);
IdxArray=0;

%% Calculation
if length(GrSize) == 1
        GrLon=Lon(1,1): GrSize(1):Lon(1,2);
        GrLat=Lat(1,2):-GrSize(1):Lat(1,1);
        NumGr=length(GrLon)*length(GrLat);
        [LatTmp,LonTmp]=meshgrid(GrLat,GrLon);
        Grid(1:NumGr,1)=reshape(LonTmp,NumGr,1);
        Grid(1:NumGr,2)=reshape(LatTmp,NumGr,1);        
else
    for i=1:length(GrSize)
        GrLon=Lon(i,1): GrSize(i):Lon(i,2);
        GrLat=Lat(i,2):-GrSize(i):Lat(i,1);
        NumGr=length(GrLon)*length(GrLat);
        GrTmp=zeros(NumGr,2);
        [LatTmp,LonTmp]=meshgrid(GrLat,GrLon);
        GrTmp(:,1)=reshape(LonTmp,NumGr,1);
        GrTmp(:,2)=reshape(LatTmp,NumGr,1);
        if i==1
            IdxRmv=[];
        elseif i==length(GrSize)
            IdxRmv=((GrTmp(:,1)>Lon(i-1,1)-Eps&GrTmp(:,1)<Lon(i-1,2)+Eps)&...
                    (GrTmp(:,2)>Lat(i-1,1)-Eps&GrTmp(:,2)<Lat(i-1,2)+Eps));
            GrTmp(IdxRmv,:)=[];
        end
        IdxArray=max(IdxArray)+1:max(IdxArray)+NumGr-sum(IdxRmv);
        Grid(IdxArray,1)=GrTmp(:,1);
        Grid(IdxArray,2)=GrTmp(:,2);
    end    
end
IdxRmv=(Grid(:,1)==0);
Grid(IdxRmv,:)=[];

end

