ti = 192;
xi = 100;

pcolor(1:301, zm, squeeze(Tf(xi,:,:,ti)).'); shading interp
hold on
contour(1:301,zm, squeeze(Tf(xi,:,:,ti)).', [17 19],'k');
hold off

%%
ti = 10;
zi = 30;
cl = [-1 1].*5e-5;

subplot(1,5,1)
pcolor(squeeze(DvDt(:,:,zi, ti)).'); shading interp
colorbar;
set(gca, 'clim', cl);

subplot(1,5,2)
pcolor(squeeze(VADV(:,:,zi, ti)).'); shading interp
colorbar;
set(gca, 'clim', cl);

subplot(1,5,3)
pcolor(squeeze(VCOR(:,:,zi, ti)).'); shading interp
colorbar;
set(gca, 'clim', cl);

subplot(1,5,4)
pcolor(-squeeze(VPRE(:,:,zi, ti)).'); shading interp
colorbar;
set(gca, 'clim', cl);

subplot(1,5,5)
tot = DvDt + VADV + VCOR + VPRE;
pcolor(squeeze(tot(:,:,zi, ti)).'); shading interp
colorbar;
set(gca, 'clim', cl);

%%
ti = 10;
zi = 20;
cl = [-1 1].*1e-8;

subplot(1,3,1)
pcolor(squeeze(DbDt(:,:,zi,ti)).'); shading interp
colorbar
set(gca, 'clim', cl);

subplot(1,3,2)
pcolor(squeeze(BADV(:,:,zi,ti)).'); shading interp
colorbar
set(gca, 'clim', cl);

subplot(1,3,3)
tot = DbDt + BADV;
pcolor(squeeze(tot(:,:,zi,ti)).'); shading interp
colorbar
set(gca, 'clim', cl);