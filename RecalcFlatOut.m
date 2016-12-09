tic
currentDirectory = pwd;
[upperPath, deepestFolder, ~] = fileparts(currentDirectory) ;
IDString=deepestFolder;
load([IDString, '_OutputsFull.mat']);

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

THETA = outputFull.T;
Q0 = outputFull.Q;
slice={0, 0, 0, tslice};
sliceEta={0,0,[1 1],tslice};
time = T(tslice(1):tslice(2))./86400; %in days
tind = floor((tslice(2)-tslice(1))/2);
 nt = length(time);
 
GenerateTheoryScalings;


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
toc