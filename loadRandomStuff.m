WTH = -TtoB.*GetVar(statefile, 'budg.nc', {'WTHMASS', '(1)'}, slice)./dh(2);;

%%
KPPGTH = GetVar(statefile, 'kppdiags.nc', {'KPPg_TH', '(1)'}, slice)./dh(2);

%%
KPPGTHz = Drv(repmat(metric, [1 1 1 nt]), KPPGTH, 'z');

%%
STRAIN = GetVar(statefile, diagfile, {'VVEL', 'UVEL', 'Dx(1)+Dy(2)'}, slice);

%%
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
%%
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

%%
DIV = DPeriodic(V, dy, 'y');

%%
WB = GetVar(statefile, 'budg.nc', {'ADVr_TH', 'Dz(1)'}, slice);

%%
bz = TtoB.*Drv(repmat(metric, [1 1 1 nt]), THETA, 'z');

%%
Vz = Drv(repmat(metric, [1 1 1 nt]), V, 'z');
by = DPeriodic(THETA, dy, 'y');

%%
ZETA = GetVar(statefile, extrafile, {'momVort3'});