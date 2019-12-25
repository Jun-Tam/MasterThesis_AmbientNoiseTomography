function OutSeismBase(Grid,Layer,VsRef,VsStrMdl,PathOut)

Depth=[0;cumsum(Layer(1:end-1))];
NumNode=size(Grid,1);

Fid=fopen(PathOut,'w');
for i=1:NumNode
    Flag=all(diff(VsStrMdl(i,1:end))>0);
    if max(VsStrMdl(i,:)) < VsRef
        if Flag==1
            DepSeisBase=interp1(VsStrMdl(i,1:end),Depth(1:end),VsRef,'linear');
        else
            DepSeisBase=interp1(VsStrMdl(i,1:5),Depth(1:end),VsRef,'linear');
        end
    else
        if Flag==1
            DepSeisBase=interp1(VsStrMdl(i,1:end),Depth(1:end),VsRef,'linear');            
        else
            DepSeisBase=interp1(VsStrMdl(i,1:9),Depth(1:9),VsRef,'linear');
        end
    end
    fprintf(Fid,'%f %f %f\n',Grid(i,1:2),1e3*DepSeisBase);
end
fclose(Fid);

end

