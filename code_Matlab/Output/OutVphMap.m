function OutVphMap(VphRef,Vmap,Vrsl,Freq,Grid,DmpTmg,Fit,Cmp,PathOutput)


%% ------------------------------------------------------------------------
if exist(PathOutput,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOutput,'s');
end
mkdir(PathOutput);

%% Velocity map
NumFreq=size(Vmap,2);
NumNode=size(Grid,1);
for i=1:NumFreq
    VphPtb=(Vmap(:,i)/VphRef(i)-1)*100;
    Fid(1)=fopen(strcat(PathOutput,num2str(Freq(i),'%2.3f'),'Hz.txt'),'w');
    for j=1:NumNode
        fprintf(Fid(1),'%.2f %.2f %.3f %.3f %.3f\n',...
            Grid(j,1:2),Vmap(j,i),Vrsl(j,i),VphPtb(j));
    end
end

%% Damping Parameters
if isempty(DmpTmg)==0
    NumDmp=length(DmpTmg);    
    for i=1:NumFreq
        OptIdx=FindOptPt(Cmp(:,i),Fit(:,i),DmpTmg);
        DmpOpt=DmpTmg(OptIdx);        
        Fid(2)=fopen(strcat(PathOutput,num2str(Freq(i),'%2.3f'),'Hz_Dmp.txt'),'w');
        fprintf(Fid(2),'%.1f %f %f\n',DmpOpt,Fit(OptIdx,i),Cmp(OptIdx,i));
        for j=1:NumDmp
            fprintf(Fid(2),'%.1f %f %f\n',DmpTmg(j),Fit(j,i),Cmp(j,i));
        end
    end
    fclose('all');
end


%% PathDensity
% Fid(2)=fopen(strcat(FileName,'_Cvr.txt'),'w');
% for j=1:NumNode
%     fprintf(Fid(2),[repmat('%.3f ',1,3),'\n'],Grid(j,1:2),PthDnsCvr(j,i));
% end

end
