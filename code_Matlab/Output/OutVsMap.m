function OutVsMap(VsStr0,Vmap,VsRsl,VphRslMdl,Layer,Grid,PathOutput)

Depth=[0;cumsum(Layer(1:end-1))];
NumLyr=length(Layer);
NumNode=size(Grid,1);

VphRsl=repmat(VphRslMdl,1,NumLyr);

%% ------------------------------------------------------------------------
if exist(PathOutput,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOutput,'s');
end
mkdir(PathOutput);

%% Vs data
for i=1:NumLyr
    Fid=fopen(strcat(PathOutput,num2str(Depth(i),'%0.1f'),'km.txt'),'w');
    Vptb=(Vmap(:,i)/VsStr0(i)-1)*100;
    for j=1:NumNode
        fprintf(Fid,['%.2f %.2f ',repmat('%.3f ',1,4),'\n'],...
            Grid(j,:),Vmap(j,i),Vptb(j),VsRsl(j,i),VphRsl(j,i));
    end
    fclose(Fid);
end

% Aftershocks
% % 3 month after the main shock (Determined by Okada et al., 2012)
% Fid=fopen('hypo.txt','r');
% Dataset=textscan(Fid,repmat('%f ',1,13),'HeaderLines',1);
% fclose(Fid);
% LatPS=Dataset{7};
% LonPS=Dataset{8};
% DepPS=Dataset{9};
% MagPS=Dataset{10};
% NumPost=length(LatPS);

% 2008-2012 (Determined by Okada et al., 2012)
% Fid=fopen('hypo.all.txt','r');
% Dataset=textscan(Fid,['%s',repmat('%f ',1,14)]);
% fclose(Fid);
% LatPS=Dataset{8};
% LonPS=Dataset{9};
% DepPS=Dataset{10};
% MagPS=Dataset{11};
% NumPost=length(LatPS);

% 2008-2010 (Determined by Ambiet noise tomography)
Fid=fopen('tomoDD.reloc.2008-2010.txt','r');
Dataset=textscan(Fid,repmat('%f ',1,24));
fclose(Fid);
LatPS1=Dataset{2};
LonPS1=Dataset{3};
DepPS1=Dataset{4};
MagPS1=Dataset{17};    % Minimum magnitude 1.5
NumPost1=length(LatPS1);

% % 2011-2012 (Determined by Ambiet noise tomography)
Fid=fopen('tomoDD.reloc.2011-2012.txt','r');
Dataset=textscan(Fid,repmat('%f ',1,24));
fclose(Fid);
LatPS2=Dataset{2};
LonPS2=Dataset{3};
DepPS2=Dataset{4};
MagPS2=Dataset{17};    % Minimum magnitude 1.5
NumPost2=length(LatPS2);
for i=1:NumLyr
    Fid(1)=fopen(strcat(PathOutput,num2str(Depth(i),'%0.1f'),'km_post1.txt'),'w');
    Fid(2)=fopen(strcat(PathOutput,num2str(Depth(i),'%0.1f'),'km_post2.txt'),'w');
    for j=1:NumPost1
        if  LatPS1(j) > min(Grid(:,2)) && LatPS1(j) < max(Grid(:,2)) && ...
            LonPS1(j) > min(Grid(:,1)) && LonPS1(j) < max(Grid(:,1)) && ...
            abs(DepPS1(j)-Depth(i)) < 0.25 && MagPS1(j) >= 1.5
            fprintf(Fid(1),[repmat('%0.3f ',1,3),'%0.1f\n'],...
                LonPS1(j),LatPS1(j),DepPS1(j),MagPS1(j));
        end
    end
    for j=1:NumPost2
        if  LatPS2(j) > min(Grid(:,2)) && LatPS2(j) < max(Grid(:,2)) && ...
            LonPS2(j) > min(Grid(:,1)) && LonPS2(j) < max(Grid(:,1)) && ...
            abs(DepPS2(j)-Depth(i)) < 0.25 && MagPS2(j) >= 0.5
            fprintf(Fid(2),[repmat('%0.3f ',1,3),'%0.1f\n'],...
                LonPS2(j),LatPS2(j),DepPS2(j),MagPS2(j));
        end
    end    
    fclose('all');
end

%% Reflector
Fid=fopen('RefInfo.txt','r');
Info=textscan(Fid,repmat('%f ',1,9));
fclose(Fid);
Data=cell2mat(Info);

for i=1:NumLyr
    Fid=fopen(strcat(PathOutput,num2str(Depth(i),'%0.1f'),'km_ref.txt'),'w');    
    for j=1:length(Data)
        if abs(Data(j,5)-Depth(i)) < 0.5
           fprintf(Fid,[repmat('%.3f ',1,5),'\n'],Data(j,[1,3,5,7,8]));
        end
    end
    fclose('all');
end

end
