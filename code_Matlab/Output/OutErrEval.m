function OutErrEval(RdErrStd,Freq,PathOut)

%% Make Directory
if exist(PathOut,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOut,'s');
end
mkdir(PathOut);

%% Output
NumFreq=length(Freq);
for i=1:NumFreq
    Fid=fopen(strcat(PathOut,num2str(Freq(i),'%2.3f'),'Hz.txt'),'w');
    IdxNonZero=RdErrStd(:,i)~=0;
    for j=1:sum(IdxNonZero)
        fprintf(Fid,'%f \n',RdErrStd(j,i));
    end
    fclose(Fid);
end

end
