% regenerateTheoryScalings
clear all; close all; clc;
addpath('~/PV_MATLAB');
addpath('~/PVINJECT_MATLAB/HelperFiles');

currentDirectory = pwd;
[upperPath, deepestFolder, ~] = fileparts(currentDirectory) ;
IDString=deepestFolder;
saveflag = 0;
statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc'; extrafile = 'extra.nc';
kppfile = 'kppdiags.nc';

TtoB = 9.81.*2e-4;
f0 = 1e-4;

Q0 = ncread(etanfile, 'TFLUX');
X = ncread(statefile, 'X');
Y = ncread(statefile, 'Y');
Z = ncread(statefile, 'Z');
Zl = ncread(statefile, 'Zl');
T = ncread(diagfile, 'T');

nx = length(X); ny = length(Y); nz = length(Z);
dx = X(2)-X(1)
dy = Y(2)-Y(1)
dz = Z(1)-Z(2) %surface only, XX-should track this through the code to ensure correct.
ts = T(2)-T(1)
% ts = 3600
dh = diff([Zl; -300]);

tslice = [1 length(T)-1];
nt = tslice(end);
gridvol = permute(repmat( dx.*dy.*abs(dh), [1 nx ny nt]), [2 3 1 4]);

    

slice={0, 0, 0, tslice};
sliceEta={0,0,[1 1],tslice};
sliceSm={0, 0, [1 6], tslice};
time = T(tslice(1):tslice(2))./86400; %in days
tind = floor((tslice(2)-tslice(1))/2);
 nt = length(time);
 
THETA = GetVar(statefile, diagfile, {'THETA', '(1)'}, sliceSm);
 GenerateTheoryScalings
 
 
% FLAT VARIABLES

output2.dJbdavg = Jbdavg;
output2.dJfdavg = Jfdavg;
output2.time = time;

FigString = [IDString, '_OutputsFlat_2.mat'];
save(FigString, 'output2', '-v7.3');
% end