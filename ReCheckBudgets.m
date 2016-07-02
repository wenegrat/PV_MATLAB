%%
%Params
 statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
% Parameters            
% dx = 1000; dy = dx; dz = 2.5;
divstrh='1e6'; divstrz='2500';
% nx = 96; ny=192; nz=200;
% ts = 3600;
tslice = [10 110];
slice={0, 0, 0, tslice};%100 120
sliceEta={0,0,[1 1],tslice};%251 271
%%
% Check if UDIAG1 equals flux terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UDIAG = GetVar(statefile, diagfile, {'UDIAG1', '(1)'}, slice);
ADVx_TH = GetVar(statefile, diagfile, {'ADVx_TH', ['Dx(1)/',divstrz]}, slice);
ADVy_TH = GetVar(statefile, diagfile, {'ADVy_TH', ['Dy(1)/',divstrz]}, slice);
ADVr_TH = GetVar(statefile, diagfile, {'ADVr_TH', ['Dz(1)/',divstrh]}, slice);
ADVT = ADVx_TH+ADVy_TH+ADVr_TH;
%%
% Confirm budget is fully closed for buoyancy
TEND = GetVar(statefile, diagfile, {'TOTTTEND', '(1)/86400'}, slice);
%%
% Alternate advection term

U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

[nx ny nz nt ]=size(W);
ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);

THETA = GetVar(statefile, diagfile, {'THETA', '(1)'}, slice);
ADVDIR = U.*Drv(dx, THETA, 'x') + V.*Drv(dy, THETA, 'y') + W.*Drv(metric, THETA, 'z');

%%
D1 = TEND - UDIAG;
D2 = TEND + ADVDIR; %oppose signs from definitions
Jbz1 = OMEGAZ.*D1;
Jbz2 = OMEGAZ.*D2;
Jbz1 = cumtrapz(t, Jbz1, 4);
Jbz2 = cumtrapz(t, Jbz2, 4);

plot(squeeze(nanmean(nanmean(Jbz1(:,:,3,:)))));
hold on
plot(squeeze(nanmean(nanmean(Jbz2(:,:,3,:)))));
hold off
%%
j=40; k=20; d =2;
plot(squeeze(UDIAG(j,k, d,:)))
hold on
plot(squeeze(-ADVT(j,k, d,:)), '--')
plot(squeeze(TEND(j,k, d,:)))
plot(squeeze(-ADVDIR(j,k,d,:)));
hold off
legend('UDIAG', 'ADVT', 'TEND');

%%
% Confirm Mom budget closes

%%
% Decompose Mom Advective Terms
UADV = GetVar(statefile, diagfile, {'Um_Advec', '(1)'}, slice);
ADVx_Um = GetVar(statefile, diagfile, {'ADVx_Um', ['Dx(1)/',divstrz]}, slice);
ADVy_Um = GetVar(statefile, diagfile, {'ADVy_Um', ['Dy(1)/',divstrz]}, slice);
ADVr_Um = GetVar(statefile, diagfile, {'ADVrE_Um', ['Dz(1)/',divstrh]}, slice);
ADVT = -ADVx_Um-ADVy_Um-ADVr_Um;
UCORI = -GetVar(statefile, diagfile, {'Um_Cori', '(1)'},slice);
RESCOR = UADV-ADVT;
%%
UADVTERM = GetVar(statefile, diagfile, {'ADVx_Um', 'ADVy_Um', 'ADVrE_Um', 'Um_Cori', ['-Dx(1)/',divstrz, '-Dy(2)/', divstrz,'-Dz(3)/', divstrh,'+(4)']}, slice);
%%
j=40; k=48; d =10;
plot(squeeze(UADVTERM(j,k, d,:)))
hold on
plot(squeeze(UADV(j,k, d,:)), '--')
hold off