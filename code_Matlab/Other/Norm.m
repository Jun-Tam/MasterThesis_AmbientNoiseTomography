function DataNrm = Norm(Data,~)

if nargin==1
    if size(Data,1)==1 || size(Data,2)==1 % 1 dimension
        DataNrm=(Data-min(Data))/(max(Data)-min(Data));
    else % 2 dimensions
        DataNrm=(Data-min(min(Data)))/(max(max(Data))-min(min(Data)));
    end
else
    DataNrm=zeros(size(Data,1),size(Data,2));
    for i=1:size(Data,1)
        DataNrm(i,:)=(Data(i,:)-min(Data(i,:)))/(max(Data(i,:))-min(Data(i,:)));
    end
end
end

