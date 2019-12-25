function DataPrtb = Prtb(Data,val)

switch nargin
    case 1
        if size(Data,1)==1 || size(Data,2)==1 % 1 dimension
            DataPrtb=(Data/(mean(Data))-1)*100;
        else % 2 dimensions
            DataPrtb=(Data/mean(mean(Data))-1)*100;
        end
    case 2
        DataPrtb=zeros(size(Data,1),size(Data,val));
        for i=1:size(Data,val)
            DataPrtb(:,i)=(Data(:,i)/(mean(Data(:,i)))-1)*100;
        end
end

end

