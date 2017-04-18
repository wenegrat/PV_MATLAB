%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN SCRIPT - PV COMPENSATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Currently documenting code called to create AGU presentation
% Written 12/19/2016

startup.m % For path and formatting defaults

f0 = 1e-4;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% DATA PROCESSING ROUTINES
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that plots assume each new model run has had 
% saveOutputs.m run on it already.
% Run by using:
% matlab -nojvm < ~/PV_MATLAB/saveOutputs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveOutputs.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verification Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentDirectory = pwd;
[upperPath, deepestFolder, ~] = fileparts(currentDirectory) ;
IDString=deepestFolder;
load([IDString, '_OutputsFlat.mat']);
BasicVerificationPlots.m

%%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% MS FIGURE ROUTINES
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: Basic Frontal Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/home/jacob13/MITgcm/projects/PV_INJECTPROJ/input/gendata.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2: 4 Panel Volume Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First pick which pre-processed run to use:
path = '/data/thomas/jacob13/PARAMSPACE/GS_025_2F/';
load([path 'GS_025_2F_OutputsFull.mat']);
U = GetVar([path 'state.nc'], [path 'diag.nc'], {'UVEL', '(1)'}, {0, 0, [1 2], 0});
% Then Make Plots
VolumePlot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3: Q Plot Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([path 'GS_025_2F_OutputsFlat.mat']);
ZETA = GetVar([path 'state.nc'], [path 'extra.nc'], {'momVort3', '(1)'}, {0, 0, [1 1], 0});
    
QPlotSurfFields.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGO Climatology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MakeDescriptiveModeWaterPlots.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Panel Volume Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VolumePlot.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Panel Volume Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_MagB_Plot.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MLI Diabatic Fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = GetVar([path 'state.nc'], [path 'diag.nc'], {'UVEL', '(1)'}, {0, 0, 0, 0});
V = GetVar([path 'state.nc'], [path 'diag.nc'], {'VVEL', '(1)'}, {0, 0, 0, 0});
W = GetVar([path 'state.nc'], [path 'diag.nc'], {'WVEL', '(1)'}, {0, 0, 0, 0});
MLIDiabaticPlot.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theory Vs Model Scatter Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CompareTheoryModels.m % Main plot
CompareTheoryModelsDomAvg.m %Domain averaged plot
CompTheoryModelsLegend.m % Made the legend separately because of export issues


path = '/data/thomas/jacob13/PARAMSPACE/GS_025_2F/';
% ZETA = GetVar([path 'state.nc'], [path 'extra.nc'], {'momVort3', '(1)'}, {0, 0, [1 1], 0});
load([path 'GS_025_2F_OutputsFull.mat']);
% U = GetVar([path 'state.nc'], [path 'diag.nc'], {'UVEL', '(1)'}, {0, 0, [1 2], 0});
load([path 'GS_025_2F_OutputsFlat.mat']);
LapeyreAnalysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theoretical J_D_t solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AnalyticBtSolution.m
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spindown PV Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basepath = '
load('/data/thomas/jacob13/PARAMSPACE/GS_SPIN_4F_SWR/GS_SPIN_4F_SWR_OutputsFull.mat');
THETA = outputFull.T;
Zl = ncread('/data/thomas/jacob13/PARAMSPACE/GS_SPIN_4F_SWR/state.nc', 'Zl');
dx = 500; dy = 500; %Uniform grid
Qo = ncread('/data/thomas/jacob13/PARAMSPACE/GS_SPIN_4F_SWR/etan.nc', 'oceQnet');
Qsw = ncread('/data/thomas/jacob13/PARAMSPACE/GS_SPIN_4F_SWR/etan.nc', 'oceQsw');;
H = ncread('/data/thomas/jacob13/PARAMSPACE/GS_SPIN_4F_SWR/etan.nc', 'KPPhbl');
QBudgetPlotSPINDOWN

%%


%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Additional analysis (not in MS)
%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% To get the non-dimensional parameters run this:
scaleEquivBuoyancyFlux.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eddy Diabatic Flux Plot (slice + APE/KE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DiabaticFluxTermsWithKE_APE.m %May need to call this first.
ComboDiabaticEddyPlot.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROMS ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RomsAnalysis.m