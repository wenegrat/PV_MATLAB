slice = {0, 0, 0, [1 359]};
statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc'; extrafile = 'extra.nc';
kppfile = 'kppdiags.nc';
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

time = T(tslice(1):tslice(2))./86400; %in days
tind = floor((tslice(2)-tslice(1))/2);
 nt = length(time);
metric = permute(repmat(Z, [1, nx, ny, 1]), [2 3 1 4]);

%% Terms in Buoyancy Budget
b = GetVar(statefile, diagfile, {'b', '(1)'}, slice);
U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

bx = DPeriodic(b, dx, 'x');
by = DPeriodic(b, dy, 'y');
%% TURBULENCE TERMS
budgfile = 'budg.nc';
divstrh = num2str(dx*dy);
divstrz = num2str(dx);

DIFF1 = GetVar(statefile, kppfile, {'KPPg_TH', ['-Dz(1)/',divstrh]}, slice);
% DIFF = DIFF+GetVar(statefile, extrafile, {'DFrI_TH', ['-Dz(1)/',divstrh]}, slice);
disp('1')
DIFF2 = DIFF1+ GetVar(statefile, budgfile, {'DFrI_TH', ['-Dz(1)/', divstrh]}, slice);
disp('2')

DIFF3 = DIFF2+ GetVar(statefile, extrafile, {'DFxE_TH', 'DFyE_TH', ['-Dx(1)/',divstrz,'-Dy(2)/',divstrz]}, slice)./permute(repmat([diff(Zl); 300+Zl(end)], [1 nx ny nt]), [2 3 1 4]);
disp('3')
TFLUX = GetVar(statefile, etanfile, {'TFLUX', '(1)'}, {0, 0, [1 1], slice{4}});
% [nx ny nd] = size(DIFF);
%%
TFLUXF = zeros(nx, ny, nz, nt);
TFLUXF(:,:,1, :) = TFLUX(:,:,:,:);
Cw = 3994;		  
H = dz;
TFLUXF = TFLUXF./(1035*Cw*H);

DIFF = DIFF3+TFLUXF;

%%
[nx ny nz nt] = size(b);
bz = Drv( repmat(metric(:,:,1:nz), [1 1 1 nt]), b,'z');
bt = NaN(size(b));
bt(:,:,:,2:end) = diff(b, 1, 4)./ts;

ADV = U.*bx + V.*by + W.*bz;

bta = squeeze(nanmean(nanmean(nanmean(bt(:,:,2:3,:)))));
ADVa = squeeze(nanmean(nanmean(nanmean(ADV(:,:,2:3,:)))));
DIFFa = squeeze(nanmean(nanmean(nanmean(DIFF(:,:,2:3,:)))));
%%
plot(time(1:nt), bta);
hold on
plot(time(1:nt), ADVa);
plot(time(1:nt), bta + ADVa);
plot(time(1:nt), TtoB.*DIFFa)
% plot(time(1:nt), -output.dJba_t(1:nt)./(f0*nx*ny*dx*dy));
hold off

%%
ix = 50; iy = 30; iz = 2;
plot(time(1:nt), squeeze(bt(ix, iy, iz, :)));
hold on
plot(time(1:nt), squeeze(ADV(ix, iy, iz, :)));
plot(time(1:nt), squeeze(bt(ix, iy, iz, :)+ADV(ix, iy, iz,:)));
plot(time(1:nt), squeeze(TtoB.*DIFF(ix, iy, iz, :)));
hold off
%%
btaz = squeeze(nanmean(nanmean(nanmean(bt(:,:,:,150:200),4))));
ADVz = squeeze(nanmean(nanmean(nanmean(ADV(:,:,:,150:200),4))));
DIFFz = TtoB.*squeeze(nanmean(nanmean(nanmean(DIFF(:,:,:, 150:200), 4))));
plot(btaz, Z); 
hold on; 
plot(ADVz, Z); 
plot(btaz+ADVz, Z);
plot(DIFFz, Z);
hold off


%%
ti = 40;
xi = 20;
cl = [-1 1].*1e-9;
ADVv = W.*bz;j
ADVh = ADV - ADVv;
subplot(1,4,1)
% pcolor(Y./1000, Z, squeeze(bt(xi,:,:,ti)).'); shading interp
pcolor(Y./1000, Z, squeeze(ADVh(xi,:,:,ti)).'); shading interp

set(gca, 'clim', cl);

subplot(1,4,2)
pcolor(Y./1000, Z, squeeze(ADVv(xi,:,:,ti)).'); shading interp
set(gca, 'clim', cl);
subplot(1,4,3)
pcolor(Y./1000, Z, TtoB.*squeeze(DIFF(xi,:,:,ti)).'); shading interp
cl = get(gca, 'clim');
hold on
contour(Y./1000, Z, squeeze(b(xi,:,:,ti)).', 20, 'k');
quiver(Y./1000, Z, squeeze(V(xi, :,:,ti)).', 100*squeeze(W(xi,:,:,ti)).', 'k')
hold off
set(gca, 'clim', cl);

subplot(1,4,4)
pcolor(X./1000, Y./1000, TtoB.*squeeze(DIFF(:,:,2,ti)).'); shading interp
 cl = get(gca, 'clim');
hold on
yt = get(gca, 'YTick');
plot(X(xi)./1000.*ones(size(yt)), yt);
% contour(Y./1000, Z, squeeze(b(xi,:,:,ti)).', 20, 'k');
hold off
set(gca, 'clim', cl);

%% DOUBLE PLOT
ti = 40;
xi = 20;
cl = [-8 8];
conts = [-100 linspace(cl(1), cl(end), 20) 100];
Bo = 9.81*2e-4*100./(1035*3994);

norm = f0.*Bo./150;
zl = [-180 0];
xl = [35 65];
vprime = V(xi,:,:,ti); - repmat(nanmean(V(xi,:,:,ti),2), [1 ny 1 1]);
% subplot(1,2,1)
[c, h]= contourf(Y./1000, Z, -squeeze(outputFull.JBz(xi,:,:,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');

% cl = get(gca, 'clim');
hold on
 contour(Y./1000, Z, squeeze(b(xi,:,:,ti)).', 20, 'k', 'LineWidth', 1.25);
dc = 4; dcz = 2;
quiver(Y(1:dc:end)./1000, Z(1:dcz:end), squeeze(vprime(1, 1:dc:end, 1:dcz:end,1)).', 100*squeeze(W(xi,1:dc:end,1:dcz:end,ti)).', 0.5,'k', 'LineWidth', 1.5)
hold off
cb= colorbar;
set(get(cb, 'ylabel'), 'string', '$\frac{{-J_D}}{f|B_o|/H}$', 'Interpreter', 'Latex', 'Rotation', 0, 'FontSize', 22);

set(gca, 'clim', cl);
set(gca, 'ylim', zl, 'xlim', xl)
colormap(cptcmap('cool-warm.cpt'));
xlabel('y (km)');
ylabel('z (m)');
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'w');
% subplot(1,2,2)
% pcolor(X./1000, Y./1000, TtoB.*squeeze(DIFF(:,:,2,ti)).'); shading interp
% % cl = get(gca, 'clim');
% hold on
% yt = get(gca, 'YTick');
% plot(X(xi)./1000.*ones(size(yt)), yt);
% % contour(Y./1000, Z, squeeze(b(xi,:,:,ti)).', 20, 'k');
% hold off
% set(gca, 'clim', cl);
% 
