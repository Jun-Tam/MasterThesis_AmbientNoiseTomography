%% Program for Ambient Noise Tomography (ANT)
% Jun TAMURA
% 2015/01/19

%% What this script computes
% This program calculates a S-wave velocity sturcture from stacked CCFs.
% First, "SPAC.m" measures phase velocity for every station pair.
% Second, "Tomography.m" conducts surface-wave tomography to convert
% path-averaged phase velocity into phase velocity maps at each central
% frequency. Finally, "VsStructure.m" performs linearized inversion for 1-D
% Vs profile at every grid, thereby constructing pseudo- 3D Vs model.

%% Reference 
% Phase velocity measurement    Nagaoka et al. (2012)
% Surface wave tomography       Barmin et al. (2001)
% Vs inversion?                  Nataf et al. (1986)

%% Description of Functions
% SPAC          Phase velocity measurement
% ErrEval       Assessment of measurement uncertainties by cluster analysis
% Tomography    Surface-wave tomography by linear inversion
% VelStructure  Construction of Vs structure by linearized inversion 
% VelStrRcv     Resolution test (CRT and RRT)
% OutSPAC       Output for spatial distribution of cross-spectra and fitting
% OutDisp       Output for phase velocity measurement
% OutVarR       Output for variance reduction in velocity measurements
% OutVphPath    Output for phase velocity on ray paths
% OutErrEval    Output for measurement uncertainties
% OutVphMap     Output for phase velocity maps
% OutVsMap      Output for Vs structure (map view)
% OutCS         Output for Vs structure (cross-sectional view)
% OutVsVph      Output for the result of Vs inversion at each grid point
% OutStr        Output for Vs structue (Vp, Vs, Rho)

%% INPUT (SPAC.m)
% Fs            Sampling frequency of raw data
% FltBand       Bandwidth of Butterworth filter
% FltPole       Number of poles for Butterworth filter
% NumSmooth     Number of smoothing times
% FreqRange     Central frequency of phase velocity
% FreqBand      Frequency band
% DataLength    Length of waveform (Original)
% WfLength      Length of waveform (Use)
% PathIn        Path for an input directory
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
%% OUTPUT (SPAC.m)
% StPos         Positions of pairs of stations (available pairs)
% StDist        Interstation distance of station pairs (available pairs)
% StAzim        Back azimuth of station pairs (available pairs)
% VphPthMdlIn   Path-averaged phase velocity
% VphPthVr      Variance reduction in detecting a phase shift
% Pos           Position of pairs of stations (all pairs)
% Dist          Interstation distance of station pairs (all pairs)
% Azim          Back azimuth of station pairs (all pairs)
% VrRef0        Variance reduction for a reference curve (undecimated)
% VphRef0       Reference phase velocity (undecimated)
% FreqVr        Frequency point of cross-spectra
% Amp           Amplitude of cross-spectra
% SpReal        Real parts of cross-spectra
% SPImag        Imaginary parts of cross-spectra
% (Wfs)         Waveforms of cross-correlation funcitons
% (Spectra)     Cross-spectra
% Freq          Central frequencies
% VphMsr        Measured phase velocity (available pairs)
% VphMsrVr      Variance reduciton (available pairs)
% VphRef        Reference phase velocity (decimated)
% VrRef         Variance reduction for a reference curve (decimated)
% Flag          Flag for available pairs
% IdxStPair     Index numbers of station pairs

%% INPUT (ErrEval.m)
% DistRng       Distance threshold for cluster analysis
% AzimRng       Azimuthal threshold for cluster analysis
% NumThr        Minimum number of data for cluster analysis

%% OUTPUT (ErrEval.m)
% RdErrStd      Standard deviation of measurement uncertainties
% RdErrMean     Mean of measurement uncertainties at each central frequency

%% INPUT (Tomography.m)
% GrSize        Grid size (degree)
% Lat           Latitude range
% Lon           Longitude range
% PrmMdl        1:Isotropic
% FlagSect      Path division (0: same pts 1: x pts per km)
% NumSect       Number of Dividing Points for each path
% NumBinAzm     Number of azimuth bin
% DmpTmg        Damping constant (Model)
% DmpTmgCbt1    Damping constant (Cbt1: 1 cell)
% DmpTmgCbt2    Damping constant (Cbt2: 4 cells)

%% OUTPUT (Tomography.m)
% "Mdl"         Model calculated with real data
% "Cbt1"        Checkerboard resolution test (1-cell anomalies)
% "Cbt2"        Checkerboard resolution test (4-cell anomalies)
% VphPth" "Out  Path-averaged phase velocity calculated from model
% VphMap" "     Phase velocity map
% VphErr" "In   A priori error of data
% VphRsl" "     Resolution matrix
% VphCmp" "     Model complexity for optimal damping constant
% VphFit" "     Model fitting for optimal damping constant
% VphDmp" "     Damping contant
% PthDnsCvr     Path density
% PthAzmCvr     Azimuthal coverage
% Grid          Grid
% GrStep        Grid step
% Cmp" "        Model complexity for each damping constant
% Fit" "        Model fitting for each damping constant

%% INPUT (VsStructure)
% Layer         Vertical grid (Depth)
% CrDist        Correlation distance
% PriMdlErr     A priori standard deviation
% VphErr0       A priori model uncertainties
% ItrStr        Number of maximum iteration
% FlagVs0       Initial model, 1:JMA
% RcvType       Type of Resolution test
% Noise         Random noise (for test)
% DmpTmgRcv     Damping constant
% FlagRcv       1: CRT, 2:RRT
% IdxDep        CRT: '24', '46', '68', RRT: 'RRT'

%% OUTPUT (VsStructure)
% "Mdl"         Model calculated with reald data
% "Rcv"         Recovery test
% VphStr" "     Phase velocity dispersion curve fitted to data
% VphErr" "     A posteriori error
% VsStr" "      Vs structure
% VsErr" "       
% VsRsl" "
% VsCmp" "
% VsFit" "      RMS residual
% VsIdxDmp" "   Index of damping (a priori standard deviation)
% VsIdxItr      Iteration number
% VsRmsItr      
% VsStr0        Initial model
% FlagDmp" "    

%% Configuration
clear all;
close all;
Path='/Users/Jun/Dropbox/ANT/';
cd(Path);
PathIn=strcat(pwd,'/Input/');
PathOut=strcat(pwd,'/Output/');
PathOut1=strcat(PathOut,'SPAC/');
PathOut2=strcat(PathOut,'Disp/');
PathOut3=strcat(PathOut,'VphVarR/');
PathOut4=strcat(PathOut,'VphPath/');
PathOut5=strcat(PathOut,'RdErr/');
PathOut6=strcat(PathOut,'MdlVph/Map/');
PathOut7=strcat(PathOut,'CbtVph/MapIn1/');
PathOut8=strcat(PathOut,'CbtVph/MapOut1/');
PathOut9=strcat(PathOut,'CbtVph/MapIn2/');
PathOut10=strcat(PathOut,'CbtVph/MapOut2/');
PathOut11=strcat(PathOut,'MdlVs/Map/');
PathOut12=strcat(PathOut,'MdlVs/Cs/');
PathOut13=strcat(PathOut,'MdlVs/Prf/');
PathOut14=strcat(PathOut,'VsANT.txt');
PathOut15=strcat(PathOut,'RcvVs/MapIn');
PathOut16=strcat(PathOut,'RcvVs/MapOut');
PathOut17=strcat(PathOut,'RcvVs/CsIn');
PathOut18=strcat(PathOut,'RcvVs/CsOut');
addpath(genpath(pwd));
PathR=strcat(pwd,'/RayleighM');
setenv('DYLD_LIBRARY_PATH','/usr/local/bin/');

% Cmp='ZZ';
% Fs=20;
% FltBand=[1/100 1];
% FltPole=3;
% FreqRange=[0.15,0.40];
% FreqStep=0.025;
% FreqBand=0.025;
% DataLength=16385;
% WfLength=8193;
% VelGrRange=[2;4];
% VelGrStep=0.01;
% PrtbMax=20/100;
% PrtbStep=0.5/100;
% PrtbRef=5/100;
% PrtbNext=7/100;
% Sigma=3;
% NumPathMax=1000;
% Min0x=2;
% Max0x=30;
% MinFindPt=3;
% NumSmooth=10;
% VrThre=0.3;
% DistRng=5;
% AzimRng=15;
% NumThr=5;
% GrSize=0.04;
% Lat=[38.54, 39.10];
% Lon=[140.36,141.08];
% FlagSect=1;
% NumSect=1;
% NumBinAzm=8;
% DmpTmg=[5:20,100];
% DmpTmgCbt1=5;
% DmpTmgCbt2=10;
% 
% %% Phase velocity Measurement
% [StPos,StDist,StAzim,VphPthMdlIn,VphPthVr,Pos,Dist,Azim,VrRef0,VphRef0,FreqVr,...
%     AmpVr,SpReal,SpImag,~,~,Freq,VphMsr,VphMsrVr,VphRef,VrRef,Flag,IdxStPair]=...
%     SPAC(PathIn,Cmp,Fs,FltBand,FltPole,FreqRange,FreqStep,FreqBand,DataLength,...
%     WfLength,VelGrRange,VelGrStep,PrtbMax,PrtbStep,PrtbRef,PrtbNext,Sigma,...
%     NumPathMax,Min0x,Max0x,MinFindPt,NumSmooth,VrThre);
% [RdErrStd,RdErrMean]=ErrEval(VphPthMdlIn,StPos,StDist,StAzim,DistRng,AzimRng,NumThr);
% OutSPAC(Dist,Azim,VrRef0,VphRef0,FreqVr,AmpVr,SpReal,SpImag,Freq,VphMsr,VelGrRange,...
%     VelGrStep,Cmp,PathOut1);
% load('SPAC.mat');
% OutDisp(VphMsr,VphRef,VrRef,VphPthMdlIn,VphMsrVr,Pos,Freq,PathOut2);
% OutVarR(Freq,VphPthVr,StDist,PathOut3);
% OutVphPath(VphRef,VphPthMdlIn,StPos,Freq,PathOut4);
% OutErrEval(RdErrStd,Freq,PathOut5);

%% Tomography
% [~,VphPthMdlOut,~,VphMapMdl,VphRmsMdl,VphErrMdl,VphRslMdl,VphCmpMdl,...
%     VphFitMdl,VphDmpMdl,PthDnsCvr,PthAzmCvr,Grid,GrStep,CmpMdl,FitMdl]=...
%     Tomography(VphRef,VphPthMdlIn,StAzim,StDist,StPos,Freq,DmpTmg,...
%     RdErrMean,0,Lon,Lat,GrSize,NumSect,NumBinAzm,FlagSect,NumPathMax);
% [~,~,VphMapCbtIn1,VphMapCbtOut1,VphErrCbt1,VphRslCbt1,~,VphCmpCbt1,VphFitCbt1,...
%     ~,~,VphDmpCbt1,~,~,~,~,CmpCbt1,FitCbt1]=Tomography(VphRef,VphPthMdlIn,...
%     StDist,StAzim,StPos,Freq,DmpTmgCbt1,RdErrMean,1,Lon,Lat,GrSize,NumSect,...
%     NumBinAzm,FlagSect,NumPathMax);
% [~,~,VphMapCbtIn2,VphMapCbtOut2,VphErrCbt2,VphRslCbt2,~,VphCmpCbt2,VphFitCbt2,...
%     ~,~,VphDmpCbt2,~,~,~,~,CmpCbt2,FitCbt2]=Tomography(VphRef,VphPthMdlIn,...
%     StDist,StAzim,StPos,Freq,DmpTmgCbt2,RdErrMean,2,Lon,Lat,GrSize,NumSect,...
%     NumBinAzm,FlagSect,NumPathMax);
% OutVphMap(VphRef,VphMapMdl,VphRslMdl,Freq,Grid,DmpTmg,FitMdl,CmpMdl,PathOut6);
% OutVphMap(VphRef,VphMapCbtIn1,VphRslCbt1,Freq,Grid,[],[],[],PathOut7);
% OutVphMap(VphRef,VphMapCbtOut1,VphRslCbt1,Freq,Grid,[],[],[],PathOut8);
% OutVphMap(VphRef,VphMapCbtIn2,VphRslCbt2,Freq,Grid,[],[],[],PathOut9);
% OutVphMap(VphRef,VphMapCbtOut2,VphRslCbt2,Freq,Grid,[],[],[],PathOut10);

%% Vs Inversions
load('TMG.mat');
% Layer=[ones(8,1);2;];
% CrDist=Layer;
% PriMdlErr=[1,.5:-.1:.3,.25:-.025:.1];
% VphErr0=repmat(0.1,length(Freq),1);
% ItrStr=20;
% FlagVs0=2;
% [VphStrMdl,VphErrMdl,VsStrMdl,VsErrMdl,VsRslMdl,VsCmpMdl,VsFitMdl,VsIdxDmpMdl,...
%     VsIdxItr,VsRmsItr,VsStr0,FlagDmpMdl]=VelStructure(VphMapMdl,VphErr0,PthDnsCvr,...
%     Layer,Grid,GrStep,Freq,ItrStr,CrDist,PriMdlErr,VphRef,FlagVs0,PathR,[],0);
% OutVsMap(VsStr0,VsStrMdl,VsRslMdl,VphRslMdl(:,11),Layer,Grid,PathOut11);
% OutCS(VsStrMdl,VsStr0,VphRslMdl(:,11),Grid,GrStep,0.01,Layer,PathOut12);
% save('STR.mat');
% OutVsVph(Grid,GrStep,Layer,Freq,VsStrMdl,VsErrMdl,VphMapMdl,VphStrMdl,VphErrMdl,...
%     VsStr0,VphRef,VsIdxDmpMdl,VsCmpMdl,VsFitMdl,PriMdlErr,FlagDmpMdl,ItrStr,...
%     VsIdxItr,VsRmsItr,PathOut13);
% OutStr(Grid,Layer,VsStrMdl,PathOut14);

%% Resolution Test
% for RcvType=4:4
%     load('STR.mat');
%     Noise=0.02;                                   % Random noise (s)
%     switch RcvType
%         case 1; DmpTmgRcv=5;  FlagRcv=1; IdxDep='24';  % SRT (Shallow)
%         case 2; DmpTmgRcv=5;  FlagRcv=1; IdxDep='46';  % SRT (Middle)
%         case 3; DmpTmgRcv=5;  FlagRcv=1; IdxDep='68';  % SRT (Deep)
%         case 4; DmpTmgRcv=10; FlagRcv=2; IdxDep='RRT'; % RRT
%     end
%     [VphPthRcvIn,VphMapRcvIn,VsStrRcvIn]=VelStrRcv(IdxDep,StDist,StAzim,StPos,VsStrMdl,...
%         Grid,GrSize,Lon,Lat,Layer,Freq,NumSect,NumBinAzm,FlagSect,PathR,FlagRcv,Noise);
%     [~,~,~,VphMapRcvOut,~,VphRslRcv,~,VphCmpRcv,VphFitRcv,~,~,VphDmpRcv,~,~,~,~]=...
%         Tomography(VphRef,VphPthRcvIn,StDist,StAzim,StPos,Freq,DmpTmgRcv,RdErrMean,...
%         0,Lon,Lat,GrSize,NumSect,NumBinAzm,FlagSect,NumPathMax);
%     [VphStrRcvOut,VphErrRcv,VsStrRcvOut,VsErrRcv,VsRslRcv,VsCmpRcv,VsFitRcv,...
%         VsIdxDmpRcv,VsIdxItrRcv,VsRmsItrRcv,~]=VelStructure(VphMapRcvOut,VphErr0,...
%         PthDnsCvr,Layer,Grid,GrStep,Freq,ItrStr,CrDist,PriMdlErr,VphRef,FlagVs0,...
%         PathR,[],FlagRcv);
%     OutVsMap(VsStr0,VsStrRcvIn,VsRslRcv,VphRslRcv,Layer,Grid,strcat(PathOut15,IdxDep,'/'));
%     OutVsMap(VsStr0,VsStrRcvOut,VsRslRcv,VphRslRcv,Layer,Grid,strcat(PathOut16,IdxDep,'/'));
%     OutCS(VsStrRcvIn,VsStr0,VphRslMdl(:,11),Grid,GrStep,0.01,Layer,strcat(PathOut17,IdxDep,'/'));
%     OutCS(VsStrRcvOut,VsStr0,VphRslMdl(:,11),Grid,GrStep,0.01,Layer,strcat(PathOut18,IdxDep,'/'));
%     save(strcat(IdxDep,'.mat'));
% end
