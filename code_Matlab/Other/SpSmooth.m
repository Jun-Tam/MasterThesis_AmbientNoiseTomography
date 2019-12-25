function SpSm = SpSmooth(Sp,Freq,NumTimes)

NumFreq=length(Freq);
SpSm=zeros(NumFreq,1);

if NumTimes==0
    SpSm=Sp;
else
    for j=1:NumTimes
        for i=1:NumFreq
            switch i
                case 1
                    SpSm(i)=Sp(i)/2+Sp(i+1)/2;
                case NumFreq
                    SpSm(i)=Sp(i-1)/2+Sp(i)/2;
                otherwise
                    SpSm(i)=Sp(i-1)/4+Sp(i)/2+Sp(i+1)/4;
            end
        end
        Sp=SpSm;
    end
end
end

