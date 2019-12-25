function [Dist,Azim,PairName,WfsTime,WfsFreq,Pos,IdxSt]=...
    FileInput(DataLength,WfLength,Fs,FltBand,FltPole,NumSmooth,Input,Cmp)

%% Channel Table Input
FileId=fopen('Vertical.txt');
Channel=textscan(FileId,['%s %f %f %s %s %f %f %f %s',repmat('%f ',1,9)]);
fclose(FileId);
StName=Channel{4};
StLat=Channel{14};
StLon=Channel{15};
StAlt=Channel{16};
IndexHiNet=51:57;
for i=IndexHiNet
    TmpName=char(StName(i));
    StName(i)=cellstr(TmpName(3:6));
end

%% Setting
Freq=Fs/2*linspace(-1,1,WfLength);
Cmb=nchoosek(1:numel(StName),2);
NumCmb=length(Cmb);
NumFreq=length(Freq);
Crd=zeros(NumCmb,4);
Tmp=zeros(NumCmb,3);
FlagData=ones(NumCmb,1);
NameTmp=cell(NumCmb);
DataTimeTmp=zeros(NumCmb,WfLength);
DataFreqTmp=zeros(NumCmb,NumFreq);
IndexWf=(DataLength-1)/2-(WfLength-1)/2+1:(DataLength-1)/2+(WfLength-1)/2+1;
[B,A]=butter(FltPole,FltBand/(Fs/2),'bandpass');

%% Waveforms load
for j=1:NumCmb
    if     StLon(Cmb(j,1)) > StLon(Cmb(j,2)); Flag=1;  P=1; N=2;
    elseif StLon(Cmb(j,1)) < StLon(Cmb(j,2)); Flag=-1; P=2; N=1;
    end
    St1=StName(Cmb(j,1));
    St2=StName(Cmb(j,2));    
    Crd(j,1)=StLon(Cmb(j,P));
    Crd(j,2)=StLat(Cmb(j,P));
    Crd(j,3)=StLon(Cmb(j,N));
    Crd(j,4)=StLat(Cmb(j,N));
    Alt1=StAlt(Cmb(j,P));
    Alt2=StAlt(Cmb(j,N));
    [arclen,az]=distance(Crd(j,2),Crd(j,1),Crd(j,4),Crd(j,3));
    Tmp(j,1)=sqrt((arclen*2*pi/360*earthRadius)^2+(Alt1-Alt2)^2)/1000;
    Tmp(j,2)=az-180;
    Tmp(j,3)=j;
    File=dir(char(strcat(Input,Cmp,'/',St1,'_',St2,'*')));
    if isempty(File)==0
        FileName=strcat(Input,Cmp,'/',File.name);
        FileId=fopen(FileName,'r');
        Wf=fread(FileId,'single');
        WfCut=Wf(IndexWf);
        fclose(FileId);
        WfFlt=filtfilt(B,A,WfCut);
        switch Flag
            case  1; DataTimeTmp(j,:)=WfFlt(1:1:end);
            case -1; DataTimeTmp(j,:)=WfFlt(end:-1:1);
        end
        Spectrum=fftshift(fft(filter(B,A,Wf),WfLength));
        DataFreqTmp(j,:)=SpSmooth(Spectrum,Freq,NumSmooth);
        NameTmp(j)=strcat(St1,'.',St2);
        FlagData(j)=0;        
    end
end

%% Sort by distance
Tmp=sortrows(Tmp,1);
WfsTime=zeros(NumCmb,WfLength);
WfsFreq=zeros(NumCmb,NumFreq);
IdxSt=zeros(NumCmb,2);
Pos=zeros(NumCmb,4);
PairName=cell(NumCmb,1);
Dist=Tmp(:,1);
Azim=Tmp(:,2);
for j=1:NumCmb
    IdxTmp=Tmp(j,3);
    WfsTime(j,:)=DataTimeTmp(IdxTmp,:);
    WfsFreq(j,:)=DataFreqTmp(IdxTmp,:);
    PairName(j)=NameTmp(IdxTmp);
    Pos(j,1:4)=Crd(IdxTmp,1:4);
    IdxSt(j,1:2)=Cmb(IdxTmp,1:2);
end
IdxRmv=FlagData==1;
IdxSt(IdxRmv,:)=[];
Dist(IdxRmv)=[];
Azim(IdxRmv)=[];
Pos(IdxRmv,:)=[];
WfsTime(IdxRmv,:)=[];
WfsFreq(IdxRmv,:)=[];
PairName(IdxRmv)=[];

end
