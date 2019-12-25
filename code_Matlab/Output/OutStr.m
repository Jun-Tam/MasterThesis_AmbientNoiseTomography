function OutStr(Grid,Layer,VsMdl,PathOut)

Depth=[0;cumsum(Layer(1:end-1))];
NumLyr=length(Layer);
NumNode=size(Grid,1);
C=[+0.9409;+2.0947;-0.8206;+0.2683;-0.0251;];

%% Vs data
Fid=fopen(PathOut,'w');
for i=1:NumLyr
    for j=1:NumNode
        Vs=VsMdl(j,i);
        Vp=C(1)+C(2)*Vs+C(3)*Vs^2+C(4)*Vs^3+C(5)*Vs^4;        
        VpVs=Vp/Vs;
        fprintf(Fid,'%.2f %.2f %.1f %.3f %.3f %.3f\n',Grid(j,:),Depth(i),Vp,Vs,VpVs);
    end
end
fclose(Fid);

end
