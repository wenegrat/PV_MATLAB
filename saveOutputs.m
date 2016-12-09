%%
clear all; close all; clc;
addpath('~/PV_MATLAB');
addpath('~/PVINJECT_MATLAB/HelperFiles');

savelarge = true;
%%
RunQAnalysis;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL VOlUME ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isoT = [0 100]; %Needed for plot below
mask = ones(nx, ny, nz, tslice(end)-tslice(1)+1);
vol = squeeze(sum(sum(sum(mask.*gridvol))));

%INTEGRATE Q
IntegrateQTerms;
%AreaIntegrateJTerms
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JFb, ~    ] = areaIntegrateJVecs(squeeze(JFz(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JFa = JFs-JFb;
[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JBb, ~    ] = areaIntegrateJVecs(squeeze(JBz(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JBa = JBs-JBb;

%%
GenerateTheoryScalings;

%%
% Create Output File
outputFull.Q = Q;
outputFull.JFz = JFz;
outputFull.JBz = JBz;
outputFull.JAx = JAx;
outputFull.JAy = JAy;
outputFull.JAz = JAz;
outputFull.JFzn= JFzN;
outputFull.JBzn= JBzN;
outputFull.JFzh= JFzH;
outputFull.JBzh = JBzH;
outputFull.T = THETA;
outputFull.time = time;
outputFull.X = X;
outputFull.Y = Y;
outputFull.Z = Z;
outputFull.gridvol = gridvol;
outputFull.ts = ts;
outputFull.dx = dx;
outputFull.dy = dy;
%%
if savelarge
FigString = [IDString, '_OutputsFull.mat'];
save(FigString, 'outputFull', '-v7.3');
end
%%
% FLAT VARIABLES
output.Qa = Qa;
output.Qt = Qt;
output.dJf = dJFdt;
output.dJb = dJBdt;
output.Jfa = JFa;
output.Jba = JBa;
output.dJfa_t = Jftota;
output.dJba_t = Jbtotpa;
output.dJbsa = Jbsurfa;
output.dJbea = Jbeddya;
output.dJfea = Jfeddya;
output.dJfga = Jfgeoa;
output.Tsurf = squeeze(THETA(:,:,1,:));
output.Q = Q0(1,1,1,1);
output.time = time;

FigString = [IDString, '_OutputsFlat.mat'];
save(FigString, 'output', '-v7.3');