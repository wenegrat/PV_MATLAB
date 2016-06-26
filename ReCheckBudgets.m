%%
%Params
 statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
% Parameters            
dx = 1000; dy = dx; dz = 2.5;
divstrh='1e6'; divstrz='2500';
nx = 96; ny=192; nz=200;
ts = 3600;
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
j=40; k=20; d =100;
plot(squeeze(UDIAG(j,k, d,:)))
hold on
plot(squeeze(-ADVT(j,k, d,:)), '--')
plot(squeeze(TEND(j,k, d,:)))
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
j=40; k=20; d =100;
plot(squeeze(ADVx_Um(j,k, d,:)))
hold on
plot(squeeze(UADV(j,k, d,:)), '--')
hold off