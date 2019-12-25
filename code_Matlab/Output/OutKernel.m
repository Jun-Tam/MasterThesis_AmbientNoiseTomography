function OutKernel(Depth,Freq,K)
NumLyr=length(Depth);
NumFrq=length(Freq);
PathOutput='/Users/Jun/Dropbox/ANT/Output/Kernel/';
if exist(PathOutput,'dir')~=0
    warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(PathOutput,'s');
end
mkdir(PathOutput);

for i=1:NumFrq
    Fid=fopen(strcat(PathOutput,num2str(Freq(i),'%2.3f'),'Hz.txt'),'w');
    for j=1:NumLyr
        fprintf(Fid,'%.3f %.3f\n',K(NumFrq-i+1,j),Depth(j));
    end
    fclose(Fid);
end

end

