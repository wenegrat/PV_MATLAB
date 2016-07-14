%% RUN Q ANALYSIS
clc; clear all; %close all;
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
TtoB = 9.81.*2e-4;

B0 = ncread(etanfile, 'TFLUX');
X = ncread(statefile, 'X');
Y = ncread(statefile, 'Y');
Z = ncread(statefile, 'Z');
T = ncread(statefile, 'T');

nx = length(X); ny = length(Y); nz = length(Z);
dx = X(2)-X(1)
dy = Y(2)-Y(1)
dz = Z(1)-Z(2)
ts = T(2)-T(1)
gridvol = dx.*dy.*dz;


tslice = [1 599];
slice={0, 0, 0, tslice};
sliceEta={0,0,[1 1],tslice};
time = T(tslice(1):tslice(2))./86400; %in days
tind = floor((tslice(2)-tslice(1))/2);
 nt = length(time);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Q Terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CalculateQTerms;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL VOlUME ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isoT = [0 100]; %Needed for plot below
mask = ones(nx, ny, nz, tslice(end)-tslice(1)+1);
vol = gridvol.*squeeze(sum(sum(sum(mask))));

%INTEGRATE Q
IntegrateQTerms;

%AreaIntegrateJTerms
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JFb, ~    ] = areaIntegrateJVecs(squeeze(JFz(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JFa = JFs-JFb;
[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JBb, ~    ] = areaIntegrateJVecs(squeeze(JBz(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JBa = JBs-JBb;

%PLOT (QBudget, dQdt)
titleString = ['Full Volume           Surface B_0: ', num2str(squeeze(B0(1,1,1)))];
QBudgetPlot;
dQdtPlot;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISOPYCNAL ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isoT = [16.1 16.5];
mask = (THETA(:,:,:,:)>isoT(1)) & (THETA(:,:,:,:)<isoT(2));
vol = gridvol.*squeeze(sum(sum(sum(mask))));

%INTEGRATE Q
IntegrateQTerms;

%AreaIntegrateJTerms
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol); %Surface Only by impermeability
JFa = JFs;
[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
JBa = JBs;

%PLOT (QBudget, dQdt)
titleString = ['Iso Volume           Surface B_0: ', num2str(squeeze(B0(1,1,1)))];
QBudgetPlot;
dQdtPlot;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAYER ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isoT = [1 100];
mask = zeros(nx, ny, nz, nt);
depthind = 50;
mask(:,:,1:depthind,:) = 1; %z > fixed depth

% deltaCrit = 0.01./(-2e-4*1035); % Convert to a delta T criteria.
% tprime = THETA - repmat(THETA(:,:,1,:), [1, 1, nz, 1]);
% mask = tprime>deltaCrit;
vol = gridvol.*squeeze(sum(sum(sum(mask))));
% 
%INTEGRATE Q
IntegrateQTerms;
% 
% %Volume Integrals of J divergences
% areaIntegrateJDiv;

%AreaIntegrateJTerms
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JFb, ~    ] = areaIntegrateJVecs(squeeze(JFz(:,:,depthind,:)), squeeze(mask(:,:,depthind,:)), dx*dy, ts, vol);
Fric = JFs-JFb;
[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JBb, ~    ] = areaIntegrateJVecs(squeeze(JBz(:,:,depthind,:)), squeeze(mask(:,:,depthind,:)), dx*dy, ts, vol);
Dia = JBs-JBb;
[JAs, dJAdt] = areaIntegrateJVecs(squeeze(JAz(:,:,3,:)), squeeze(mask(:,:,3,:)), dx*dy, ts, vol);
[JAb, ~    ] = areaIntegrateJVecs(squeeze(JAz(:,:,depthind,:)), squeeze(mask(:,:,depthind,:)), dx*dy, ts, vol);
Adv = JAs-JAb;

%PLOT (QBudget, dQdt)
titleString = ['Layer Analysis           Surface B_0: ', num2str(squeeze(B0(1,1,1)))];
QBudgetPlotAdv;
% dQdtPlot;

%%
RiPlot;

%%
toc./60 % in minutes