function [StPos,StDist,StAzim,VphPth,VphPthVr,Pos,Dist,Azim,VrRef0,VphRef0,...
    Freq,Amp,SpReal,SpImag,Wfs,Spectra,FreqCent,VphMsr,VphMsrVr,VphRef,VrRef,...
    Flag,IdxStPair]=SPAC(PathIn,Cmp,Fs,FltBand,FltPole,FreqRange,FreqStep,...
    FreqBand,DataLength,WfLength,VelGrRange,VelGrStep,PrtbMax,PrtbStep,...
    PrtbRef,PrtbNext,Sigma,NumPathMax,Min0x,Max0x,MinFindPt,NumSmooth,VrThre)


%% Description of functions
% "FileInput.m" loads stacked CCFs.
% "RefVelCal.m" applies spatial auto-correlation method to cross-spectra 
%               and determines a referece dispersion curve and reference
%               cross-spectra.
% "DisperCal.m" measures phase velocity for every station pair by detecting 
%               a phase shift between the reference and observed cross-spectra.

%- INPUT
% Fs            Sampling frequency of raw data
% FltBand       Bandwidth of Butterworth filter
% FltPole       Number of poles for Butterworth filter
% NumSmooth     Number of smoothing times
% FreqRange     Central frequency of phase velocity
% FreqBand      Frequency band
% DataLength    Length of waveform (Original)
% WfLength      Length of waveform (Use)
% PathIn
% Cmp           Component of CCFs
% VelGrRange    Velocity range for grid search (Ref)
% VelGrStep     Velocity step for grid search (Ref)
% PrtbMax       Maximum perturbation for grid Search (Syn)
% PrtbStep      Perturbation step for grid search (Syn)
% FreqStep      Frequency step for phase velocity
% Sigma         Outliar removal
% NumPathMax    Maximum number of ray paths
% Min0x         Minimum distance boundary (xth 0 crossing)
% Max0x         Maximum distance boundary (xth 0 crossing)
% MinFindPt     Criteria (Number of found Pts)
% PrtbRef       Perturbation
% PrtbNext      Perturbation
% VrThre        Threshold for variance reduction
%
%- OUTPUT
% StPos         Positions of pairs of stations (available pairs)
% StDist        Interstation distance of station pairs (available pairs)
% StAzim        Back azimuth of station pairs (available pairs)
% VphPth        Path-averaged phase velocity
% VphPthVr      Variance reduction in detecting a phase shift
% Pos           Position of pairs of stations (all pairs)
% Dist          Interstation distance of station pairs (all pairs)
% Azim          Back azimuth of station pairs (all pairs)
% VrRef0        Variance reduction for a reference curve (undecimated)
% VphRef0       Reference phase velocity (undecimated)
% Freq          Frequency points of cross-spectra
% Amp           Amplitude of cross-spectra 
% SpReal        Real parts of cross-spectra
% SPImag        Imaginary parts of cross-spectra
% Wfs           Waveforms of cross-correlation funcitons
% Spectra       Cross-spectra
% FreqCent      Central frequencies
% VphMsr        Measured phase velocity (available pairs)
% VphMsrVr      Variance reduciton (available pairs)
% VphRef        Reference phase velocity (decimated)
% VrRef         Variance reduction for a reference curve (decimated)
% Flag          Flag for available pairs
% IdxStPair     Index numbers of station pairs

%% Dispersion measurement
tic;
[Dist,Azim,~,Wfs,Spectra,Pos,IdxSt]=...
    FileInput(DataLength,WfLength,Fs,FltBand,FltPole,NumSmooth,PathIn,Cmp);
[VphRef0,VrVphRef0,VrRef0,Amp,Freq,Sp,SpReal,SpImag]=RefVelCal(Spectra,...
    Dist,VelGrRange,VelGrStep,FreqRange,FreqStep,Fs,WfLength,PrtbRef);
[Flag,VphMsr,VphRef,VrRef,FreqCent,VphMsrVr]=...
    DisperCal(Freq,VphRef0,Amp,VrVphRef0,Sp,Dist,FreqRange,FreqStep,...
    FreqBand,PrtbMax,PrtbStep,Sigma,Min0x,Max0x,MinFindPt,PrtbNext,VrThre);

%% Output
NumData=length(Dist);
NumFreq=length(FreqCent);
StPos=zeros(NumPathMax,4,NumFreq);
StDist=zeros(NumPathMax,NumFreq);
StAzim=zeros(NumPathMax,NumFreq);
VphPth=zeros(NumPathMax,NumFreq);
VphPthVr=zeros(NumPathMax,NumFreq);
IdxStPair=zeros(NumPathMax,2,NumFreq);

for i=1:NumFreq
    k=1;
    for j=1:NumData        
        if Flag(i,j)~=0
            StPos(k,1:4,i)=Pos(j,1:4);
            StDist(k,i)=Dist(j);
            StAzim(k,i)=Azim(j);
            VphPth(k,i)=VphMsr(i,j);
            VphPthVr(k,i)=VphMsrVr(i,j);
            IdxStPair(k,1:2,i)=IdxSt(j,1:2);
            k=k+1;
        end
    end
end
LapseTime=toc;
display(strcat('SPAC:',num2str(LapseTime),'sec'));

end
