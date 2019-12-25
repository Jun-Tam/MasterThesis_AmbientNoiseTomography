function OutVarR(Freq,VphPthVr,StDist,PathOutput)

%% Make Ouput Folder
if exist(PathOutput,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOutput,'s');
end
mkdir(PathOutput);

NumFreq=length(Freq);
for i=1:NumFreq
    FileName=strcat(PathOutput,num2str(Freq(i),'%2.3f'),'Hz');    
    Fid=fopen(strcat(FileName,'_Vr.txt'),'w');
    IdxPath=VphPthVr(:,i)~=0;
    NumPath=find(IdxPath, 1,'last' );
    for j=1:NumPath
        fprintf(Fid,'%f %f\n',StDist(j,i),VphPthVr(j,i)*100);
    end
end

% %%
% X=0:1:70;
% Y=(30:100)/100;
% NumFreq=length(Freq);
% NumX=length(X);
% NumY=length(Y);
% NumBin=NumX*NumY;
% [BinX,BinY]=meshgrid(X,Y);
% BinX=reshape(BinX',NumBin,1);
% BinY=reshape(BinY',NumBin,1);
% 
% %% Output (Variance Reduction)
% for i=1:NumFreq
%     FileName=strcat(PathOutput,num2str(Freq(i),'%2.3f'),'Hz');    
%     Fid=fopen(strcat(FileName,'_Vr.txt'),'w');
%     IdxPath=VphPthVr(:,i)~=0;
%     NumPath=find(IdxPath, 1,'last' );
%     BinCnt=zeros(NumX,NumY);
%     for j=1:NumPath
%         IdxX=find(abs(X-StDist(j,i))==min(abs(X-StDist(j,i))));
%         IdxY=find(abs(Y-VphPthVr(j,i))==min(abs(Y-VphPthVr(j,i))));
%         BinCnt(IdxX,IdxY)=BinCnt(IdxX,IdxY)+1;        
%     end
%     PrbDns=reshape(BinCnt/NumPath,NumBin,1);
%     for j=1:NumBin
%         fprintf(Fid,'%f %f %f \n',BinX(j),BinY(j),PrbDns(j));
%     end    
% end
% fclose(Fid);

end

