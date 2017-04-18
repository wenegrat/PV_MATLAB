%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized Omega RHS Terms
% Following Thomas Et Al. 2010 (Subpolar front pt II), Eq. 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc'; extrafile = 'extra.nc';
kppfile = 'kppdiags.nc';

TtoB = 9.81.*2e-4;
f0 = 1e-4;
nx = 150; ny = 200; nz = 50;
dx = 500; dy = 500;
dz = 0; %unused
ts = 7200;
nt = 200;
st = 48;
% st = 200;
slice = {0, 0, 0, [st st+nt-1]};
sliceEta = {0, 0, [1 1], [st st+nt-1]};
sizes = [nx, ny, nz, nt];
ztmp = ncread(statefile, 'Z');
ztmp = ztmp(1:nz);
Z = ztmp;
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);

%% Decompose into GEO/AGEO
b = GetVar(statefile, diagfile, {'b', '(1)'}, slice);
bx = DPeriodic(b, dx,'x');
by = DPeriodic(b, dy, 'y');
bz = Drv(metric, b,'z');

ugz = - 1./f0 .* by;
vgz =   1./f0 .* bx;

eta = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.81*(1)'},sliceEta));
dEtadx = DPeriodic(eta, dx,'x');
dEtady = DPeriodic(eta, dy, 'y');

ug  = cumtrapz(Z, ugz, 3);
ug = ug+permute(repmat(dEtady./f0, [1 1 1 nz]), [1 2 4 3]);

vg  = cumtrapz(Z, vgz, 3);
vg = vg+permute(repmat(-dEtadx./f0, [1 1 1 nz]), [1 2 4 3]);

U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

ua = U-ug;
va = V-vg;

ZETA = GetVar(statefile, extrafile, {'momVort3', '(1)'}, slice);
STRAIN = GetVar(statefile, extrafile, {'Strain', '(1)'}, slice);

%% I - GEOSTROPHIC FORCING
Ugx = DPeriodic(ug, dx, 'x');
Ugy = DPeriodic(ug, dy, 'y');
Vgx = DPeriodic(vg, dx, 'x');
Vgy = DPeriodic(vg, dy, 'y');

Fgx = 2*(Ugx.*bx + Vgx.*by);
Fgy = 2*(Vgy.*by + Ugy.*bx);

%% II - FRICTIONAL FORCES
[~,~, ~, ~, FX, ~, ~, ~, ~ , ~]= returnMomTerms(diagfile, statefile, etanfile,extrafile, sizes, slice, dx, dy,dz );
[~,~, ~, ~, FY, ~]= returnMomTermsX(diagfile, statefile, etanfile,extrafile, sizes, slice, dx, dy,dz );

dzFY = Drv(metric, FY, 'z');
% ddzFY= Drv(metric, dzFY, 'z');
dzFX = Drv(metric, FX, 'z');
% ddzFX= Drv(metric, dzFX, 'z');

Fricx = f0.*dzFY;
Fricy = f0.*(-dzFX);

%% III - BUOYANCY SOURCES/SINKS
dTdt = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND',  '(1)/86400'}, slice);
TADV = U.*bx + V.*by + W.*bz;

wprime = W - repmat(nanmean(W), [nx 1 1 1]);
bprime = b - repmat(nanmean(b), [nx 1 1 1]);
wpbp = wprime.*bprime;

Q = dTdt+ TADV;
dQdx = DPeriodic(Q, dx, 'x');
dQdy = DPeriodic(Q, dy, 'y');

Diax = -dQdx;
Diay = -dQdy;

%% IV - THERMAL WIND IMBALANCE TENDENCY
[~, ~, vaz, ~] = gradient(va, 1, 1, Z, 1);
[~, ~, uaz, ~] = gradient(ua, 1, 1, Z, 1);

[~, ~, ~, dVazdt] = gradient(vaz, ts);
[~, ~, ~, dUazdt] = gradient(uaz, ts);

Vazx = DPeriodic(vaz, dx, 'x');
Vazy = DPeriodic(vaz, dy, 'y');
Uazx = DPeriodic(uaz, dx, 'x');
Uazy = DPeriodic(uaz, dy, 'y');

TWIx = f0.*(-f0.*dVazdt + -f0.*ug.*Vazx + -f0.*vg.*Vazy);
TWIy = f0.*(f0.*dUazdt + f0.*ug.*Uazx + f0.*vg.*Uazy);

%% V - LARGE RO FLOW
vagx = DPeriodic(va, dx, 'x');
vagy = DPeriodic(va, dy, 'y');
uagx = DPeriodic(ua, dx, 'x');
uagy = DPeriodic(ua, dy, 'y');

[~, ~, dzv, ~] = gradient( ua.*vagx + va.*vagy, 1, 1, Z, 1); %Supposed to have w?
[~, ~, dzzv, ~] = gradient( dzv, 1, 1, Z, 1); %Supposed to have w?

[~, ~, dzu, ~] = gradient( ua.*uagx + va.*uagy, 1, 1, Z, 1); %Supposed to have w?
[~, ~, dzzu, ~] = gradient( dzu, 1, 1, Z, 1); %Supposed to have w?

Rx = -f0.*dzzv;% ddz( uag dot nabla vag);
Ry =  f0.*dzzu;% ddz( uag dot nabla uag);

%%
[omega, phi,theta, ub, vb] =calculateOmegaMITgcm(ug, vg, dx);


%%
ti = 5;
cl = [-1 1].*1e-10;
cv = linspace(cl(1), cl(end), 20);
figure

% Geostrophic

subplot(5,2, 1)
[c, h] = contourf(1:ny, Z, squeeze(Fgx(2,:,:,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny, Z, squeeze(b(2,:,:,ti)).', 'k');
hold off
title('I - Geo - X');

subplot(5,2, 2)
[c, h] = contourf(1:ny, Z, squeeze(Fgy(2,:,:,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny, Z, squeeze(b(2,:,:,ti)).', 'k');
hold off
title('I - Geo - Y');

% Friction
subplot(5,2, 3)
[c, h] = contourf(1:ny, Z, squeeze(Fricx(2,:,:,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny, Z, squeeze(b(2,:,:,ti)).', 'k');
hold off
title('II - Fric - X');

subplot(5,2, 4)
[c, h] = contourf(1:ny, Z, squeeze(Fricy(2,:,:,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny, Z, squeeze(b(2,:,:,ti)).', 'k');
hold off
title('II - Fric - Y');

% Buoyancy
subplot(5,2, 5)
[c, h] = contourf(1:ny, Z, squeeze(Diax(2,:,:,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny, Z, squeeze(b(2,:,:,ti)).', 'k');
hold off
title('III - Dia - X');

subplot(5,2, 6)
[c, h] = contourf(1:ny, Z, squeeze(Diay(2,:,:,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny, Z, squeeze(b(2,:,:,ti)).', 'k');
hold off
title('III - Dia - Y');

% TWI
subplot(5,2, 7)
[c, h] = contourf(1:ny, Z, squeeze(TWIx(2,:,:,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny, Z, squeeze(b(2,:,:,ti)).', 'k');
hold off
title('IV - TWI - X');

subplot(5,2, 8)
[c, h] = contourf(1:ny, Z, squeeze(TWIy(2,:,:,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny, Z, squeeze(b(2,:,:,ti)).', 'k');
hold off
title('IV - TWI - Y');

% ROSSBY
subplot(5,2, 9)
[c, h] = contourf(1:ny, Z, squeeze(Rx(2,:,:,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny, Z, squeeze(b(2,:,:,ti)).', 'k');
hold off
title('V - Rossby - X');

subplot(5,2, 10)
[c, h] = contourf(1:ny, Z, squeeze(Ry(2,:,:,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny, Z, squeeze(b(2,:,:,ti)).', 'k');
hold off
title('V - Rossby - Y');
%%
figure
subplot(2,2,1)
[c, h] = contourf(1:nx, 1:ny, squeeze(Fgx(:,:,1,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:nx, 1:ny, squeeze(b(:,:,1,ti)).', 'k');
hold off
title('I - Geo - X');

subplot(2,2,2)
[c, h] = contourf(1:nx, 1:ny, squeeze(Fgy(:,:,1,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:nx, 1:ny, squeeze(b(:,:,1,ti)).', 'k');
hold off
title('I - Geo - Y');

subplot(2,2,3)
[c, h] = contourf(1:nx, 1:ny, squeeze(Fricx(:,:,4,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:nx, 1:ny, squeeze(b(:,:,1,ti)).', 'k');
hold off
title('III - Fric - X');
subplot(2,2,4)
[c, h] = contourf(1:nx, 1:ny, squeeze(Fricy(:,:,4,ti)).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:nx, 1:ny, squeeze(b(:,:,1,ti)).', 'k');
hold off
title('III - Fric - Y');


%%
figure
subplot(2,2,1)
[c, h] = contourf(1:ny, Z, squeeze(nanmean(Fgx(:,:,:,ti))).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny,Z, squeeze(nanmean(b(:,:,:,ti))).', 'k');
hold off
title('I - Geo - X');

subplot(2,2,2)
[c, h] = contourf(1:ny, Z, squeeze(nanmean(Fgy(:,:,:,ti))).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny,Z, squeeze(nanmean(b(:,:,:,ti))).', 'k');
hold off
title('I - Geo - Y');

subplot(2,2,3)
[c, h] = contourf(1:ny, Z, squeeze(nanmean(Fricx(:,:,:,ti))).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny,Z, squeeze(nanmean(b(:,:,:,ti))).', 'k');
hold off
title('III - Fric - X');

subplot(2,2,4)
[c, h] = contourf(1:ny, Z, squeeze(nanmean(Fricy(:,:,:,ti))).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:ny,Z, squeeze(nanmean(b(:,:,:,ti))).', 'k');
hold off

title('III - Fric - Y');

%%
zetas = ZETA - omega;

cv = linspace(-2, 2, 20);
cl = [cv(1) cv(end)];
figure
subplot(2,1,1)
[c, h] = contourf(1:nx, 1:ny, squeeze(omega(:,:,4,ti)./f0).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:nx, 1:ny, squeeze(b(:,:,1,ti)).', 'k');
hold off

subplot(2,1,2)
[c, h] = contourf(1:nx, 1:ny, squeeze(zetas(:,:,4,ti)./f0).', cv);
set(h, 'edgecolor', 'none')
set(gca, 'clim', cl);
hold on
contour(1:nx, 1:ny, squeeze(b(:,:,1,ti)).', 'k');
hold off
%%
magb = sqrt(bx.^2 + by.^2);
% magb = -bx;
% magb = Fricy;
zl = 1:3;
fmag = sqrt(Fricx.^2 + Fricy.^2);
% fmag = Fricy;
geomag = sqrt(Fgy.^2 + Fgx.^2);

zvec = reshape(ZETA(:,:,zl,:)./f0, [nx*ny*length(Z(zl))*nt, 1]);
fricvec = reshape(fmag(:,:,zl,:), [nx.*ny*length(Z(zl)).*nt, 1]);
geovec = reshape(geomag(:,:,zl,:), [nx*ny*length(Z(zl)).*nt, 1]);
strainvec = reshape(STRAIN(:,:,zl,:)./f0, [nx*ny*length(Z(zl)).*nt, 1]);

bvec = reshape(f0.*magb(:,:,zl,:), [nx*ny*length(Z(zl)).*nt, 1]);

%%
figure
subplot(2,1,1)
scatter(fricvec, bvec);
mask = isfinite(bvec+fricvec);

cr = corr(bvec(mask), fricvec(mask));
title(num2str(cr))
onetoone


subplot(2,1,2);
scatter(geovec, zvec);

%%
figure
scatter(squeeze(nanmean(nanmean(nanmean(fmag(:,:,zl,:))))), 1.*squeeze(nanmean(nanmean(nanmean(f0.*magb(:,:,zl,:))))))
onetoone