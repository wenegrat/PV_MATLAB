%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISOPYCNAL ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
THETA = outputFull.T;
gridvol = outputFull.gridvol;
JFz = outputFull.JFz;
JBz = outputFull.JBz;
dx = outputFull.dx; dy = outputFull.dy;
ts = outputFull.ts;
Q = outputFull.Q;
X = outputFull.X;
Y = outputFull.Y;
Z = outputFull.Z;
time = outputFull.time;
T = time;
%%
isoT = [16.8 17.2];
mask = (THETA(:,:,:,:)>isoT(1)) & (THETA(:,:,:,:)<isoT(2));
vol = squeeze(sum(sum(sum(mask.*gridvol))));

%INTEGRATE Q
IntegrateQTerms;

%AreaIntegrateJTerms
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol); %Surface Only by impermeability
JFa = JFs;
[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
JBa = JBs;

%PLOT (QBudget, dQdt)
titleString = ['Iso Volume           Surface B_0: ', num2str(squeeze(Q0(1,1,1)))];
QBudgetPlot;
dQdtPlot;
