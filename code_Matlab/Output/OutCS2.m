function OutCS2(VlcPos1,VlcPos2,Name,CsDsp1,CsDsp2,Grid1,Grid2,Layer,Drct1PS1,...
    Drct2PS1,DepPS1,MagPS1,Drct1PS2,Drct2PS2,DepPS2,MagPS2,AltMap,VsMap,VsPtb,...
    VphRsl,RefPos1,RefPos2,RefDep,RefStr,RefDip,StPos1,StPos2,PathOutput)

%-------------------------------------------------------------------------%
NumPost1=length(MagPS1);
NumPost2=length(MagPS2);
NumVlc=length(VlcPos1);
NumRef=length(RefPos1);
NumSt=length(StPos1);
Depth=[0;cumsum(Layer(1:end-1))];
[X,Y]=meshgrid(CsDsp2,Depth);
X=reshape(X',numel(X),1);
Y=reshape(Y',numel(X),1);
Z1=reshape(VsMap,length(CsDsp1),numel(X))';
Z2=reshape(VsPtb,length(CsDsp1),numel(X))';
CsDsp2Mesh=linspace(min(CsDsp2),max(CsDsp2),10000);
%-------------------------------------------------------------------------%

mkdir(strcat(PathOutput,Name));
for j=1:length(CsDsp1)
    %% File Open
    OutDirName=strcat(PathOutput,Name,num2str(CsDsp1(j),'%.2f'));
    Fid(1)=fopen(strcat(OutDirName,'_vel.txt'),'w');
    Fid(2)=fopen(strcat(OutDirName,'_rsl.txt'),'w');    
    Fid(3)=fopen(strcat(OutDirName,'_grd.txt'),'w');    
    Fid(4)=fopen(strcat(OutDirName,'_elv.txt'),'w');
    Fid(5)=fopen(strcat(OutDirName,'_vlc.txt'),'w');
    Fid(6)=fopen(strcat(OutDirName,'_ref.txt'),'w');
    Fid(7)=fopen(strcat(OutDirName,'_St.txt'),'w');  
    Fid(8)=fopen(strcat(OutDirName,'_pst1.txt'),'w');
    Fid(9)=fopen(strcat(OutDirName,'_pst2.txt'),'w');
    
    %% Swave Velocity
    for i=1:numel(X)
        fprintf(Fid(1),'%.3f %.3f %.3f %.3f\n',X(i),Y(i),Z1(i,j),Z2(i,j));
    end
    
    %% Resolution
    RslThrsh=0.3;
    [RslMax,~]=max(VphRsl(j,:));    
    if RslMax > RslThrsh
        VphRslIntp=interp1(CsDsp2,VphRsl(j,:),CsDsp2Mesh);
        IdxCross=find(VphRslIntp - RslThrsh > 0);
        if length(IdxCross) >= 2
            RslX1=CsDsp2Mesh(min(IdxCross));
            RslX2=CsDsp2Mesh(max(IdxCross));
            if abs(RslX1-RslX2) > 0.04
                fprintf(Fid(2),'%.4f %.4f\n',RslX1,0);
                fprintf(Fid(2),'%.4f %.4f\n',RslX1,10);
                fprintf(Fid(2),'%s \n','>');
                fprintf(Fid(2),'%.4f %.4f\n',RslX2,0);
                fprintf(Fid(2),'%.4f %.4f\n',RslX2,10);
                fprintf(Fid(2),'%s \n','>');
            end
        end
    end
    
    %% Grid
    [GrX,GrY]=meshgrid(Grid2(abs(Grid1-CsDsp1(j)) < 1e-4),Depth);
    GrX=reshape(GrX,numel(GrX),1);
    GrY=reshape(GrY,numel(GrY),1);
    for i=1:length(GrX)
        fprintf(Fid(3),'%.2f %.2f\n',GrX(i),GrY(i));
    end
    
    %% Elevation
    Elevation=AltMap(:,CsDsp1(j)==CsDsp1);
    for i=1:length(Elevation)
        fprintf(Fid(4),'%.2f %.3f\n',CsDsp2(i),Elevation(i)/1000);
    end
    
    %% Volcanoes
    for i=1:NumVlc
        if abs(VlcPos2(i)-CsDsp1(j)) < 0.05
            fprintf(Fid(5),'%.3f %.1f\n',VlcPos1(i),.05);
        end
    end
    
    %% Reflectors
    for i=1:NumRef
        if abs(RefPos2(i)-CsDsp1(j)) < 0.04 && max(Depth) > RefDep(i) && ...
           RefPos1(i) > min(CsDsp2) && RefPos1(i) < max(CsDsp2)
            fprintf(Fid(6),[repmat('%.3f ',1,4),'\n'],...
                RefPos1(i),RefDep(i),RefDip(i),RefStr(i));           
        end
    end
    
    %% Stations
    for i=1:NumSt
        if abs(StPos1(i)-CsDsp1(j)) < 0.04
            fprintf(Fid(7),'%.3f %.1f\n',StPos2(i),.6);
        end
    end

    %% Postshocks1
    for i=1:NumPost1
        if abs(Drct1PS1(i)-CsDsp1(j)) < 0.01 && MagPS1(i) >= 1.5 && ...
                min(Grid1) < Drct1PS1(i) && max(Grid1) > Drct1PS1(i) &&...
                min(Grid2) < Drct2PS1(i) && max(Grid2) > Drct2PS1(i)            
            fprintf(Fid(8),'%.3f %.3f %.1f\n',Drct2PS1(i),DepPS1(i),MagPS1(i));
        end
    end

    %% Postshocks2
    for i=1:NumPost2
        if abs(Drct1PS2(i)-CsDsp1(j)) < 0.01 && MagPS2(i) >= 0.5 && ...
                min(Grid1) < Drct1PS2(i) && max(Grid1) > Drct1PS2(i) &&...
                min(Grid2) < Drct2PS2(i) && max(Grid2) > Drct2PS2(i)
            fprintf(Fid(9),'%.3f %.3f %.1f\n',Drct2PS2(i),DepPS2(i),MagPS2(i));
        end
    end    
    fclose('all');    
end

end
