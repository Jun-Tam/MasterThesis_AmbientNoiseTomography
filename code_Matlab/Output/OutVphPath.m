function OutVphPath(VphRef,VphPthMdlIn,StPos,Freq,PathOutput)

%% Make Ouput Folder
if exist(PathOutput,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOutput,'s');
end
mkdir(PathOutput);

%% Make Color Palette
PltSeism=fopen('PltSeism.cpt','r');
PltJet=fopen('PltHot.cpt','r');
RgbInfo1=cell2mat(textscan(PltSeism,repmat('%f ',1,8),'HeaderLines',3));
RgbInfo2=cell2mat(textscan(PltJet,repmat('%f ',1,8),'HeaderLines',3));
fclose('all');

NumFreq=length(Freq);
NumClr1=length(RgbInfo1(:,1));
NumClr2=length(RgbInfo2(:,1));
RngSeism=RgbInfo1(:,1);
% RngJet=RgbInfo2(:,1);
ClrPlts(1:NumClr1,1:3,1)=RgbInfo1(:,2:4)/255;
ClrPlts(1:NumClr2,1:3,2)=RgbInfo2(:,2:4)/255;
% ClrPlts(1:NumClr1,1:3,1)=interp1(0:1/NumPlt1:1-1/NumPlt1,ClrSeism(:,1:3),...
%                                  0:1/NumClr1:1-1/NumClr1,'linear','extrap');
% ClrPlts(1:NumClr2,1:3,2)=interp1(0:1/NumPlt2:1-1/NumPlt2,ClrJet(:,1:3),...
%                                  0:1/NumClr2:1-1/NumClr2,'linear','extrap');

%% Output (Path)
for i=1:NumFreq
    FileName=strcat(PathOutput,num2str(Freq(i),'%2.3f'),'Hz');
    Fid(1)=fopen(strcat(FileName,'_Vel.txt'),'w');
%     Fid(2)=fopen(strcat(FileName,'_Vr.txt'),'w');
    Fid(3)=fopen(strcat(FileName,'_Pth.txt'),'w');
    IdxPath=VphPthMdlIn(:,i)~=0;
    NumPath=find(IdxPath, 1,'last' );
    for j=NumPath:-1:1
        PthPrtb=(VphPthMdlIn(j,i)/VphRef(i)-1)*100;
        %% Seismic Palette (Perturbation)
        if PthPrtb > max(RngSeism)
            IdxClr1=NumClr1;
        elseif PthPrtb < min(RngSeism)
            IdxClr1=1;
        else
            IdxClr1=find(abs(RngSeism-PthPrtb)==min(abs(RngSeism-PthPrtb)),1);
        end
%         PthVr=VphPthVr(j,i)*100;        
%         IdxClr2=find(abs(RngJet-PthVr)==min(abs(RngJet-PthVr)),1);
%         RgbJet=round(255*ClrPlts(IdxClr2,1:3,2));
%         fprintf(Fid(2),[repmat('%.3f ',1,5),repmat('%d ',1,3),'\n'],...
%             StPos(j,:,i),PthVr,RgbJet);
        RgbSeism=round(255*ClrPlts(IdxClr1,1:3,1));
        fprintf(Fid(1),[repmat('%.3f ',1,6),repmat('%d ',1,3),'\n'],...
            StPos(j,:,i),VphPthMdlIn(j,i),PthPrtb,RgbSeism);
        fprintf(Fid(3),[repmat('%.3f ',1,2),'\n'],StPos(j,1:2,i));
        fprintf(Fid(3),[repmat('%.3f ',1,2),'\n'],StPos(j,3:4,i));
        fprintf(Fid(3),'%s \n','>');
    end
    fclose('all');
end

end

