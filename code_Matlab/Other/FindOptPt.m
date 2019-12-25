function OptInx = FindOptPt(Cmplxty,Fitting,Damping)

NumDamping=length(Damping);
TrdOff=zeros(NumDamping,3);
for i=1:NumDamping
    TrdOff(i,1)=Cmplxty(i);
    TrdOff(i,2)=Fitting(i);
    TrdOff(i,3)=i;
end
TrdOff=sortrows(TrdOff,1);
lineVec=TrdOff(end,1:2)-TrdOff(1,1:2);
lineVecN=lineVec/sqrt(sum(lineVec.^2));
vecFromFirst=bsxfun(@minus,TrdOff(:,1:2),TrdOff(1,1:2));
scalarProduct=dot(vecFromFirst,repmat(lineVecN,NumDamping,1),2);
vecFromFirstParallel=scalarProduct*lineVecN;
vecToLine=vecFromFirst-vecFromFirstParallel;
distToLine=sqrt(sum(vecToLine.^2,2));
[~,OptPt]=max(distToLine);
OptInx=TrdOff(OptPt,3);

end

