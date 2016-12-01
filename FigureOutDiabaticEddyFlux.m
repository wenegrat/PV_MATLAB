t = 100;
Jbterm = JBz - JBzN - JBzH;
Jbterm = JBz;
subplot(3,1,1)
% pcolor(squeeze(ZETA(:,:,1,t))./f0); shading interp
% pcolor(squeeze(mfk(:,:,1,t)).^2.*Hbl(:,:,1,t)./f0); shading interp
pcolor(squeeze(ZETA(:,:,2,t)./f0)); shading interp

hold on
contour(squeeze(Jbterm(:,:,2,t)),30, 'k');
hold off
colorbar;
set(gca, 'clim', [-1 1]*2);
% set(gca, 'clim', [-1 1].*1e-11)

subplot(3,1,2)
pcolor(squeeze(Jbterm(:,:,2,t))); shading interp
colorbar;
set(gca, 'clim', [-1 1].*1e-11)
hold on
contour(squeeze(Jbterm(:,:,2,t)),30, 'k');
hold off

subplot(3,1,3)
pcolor(squeeze(STRAIN(:,:,1,t))./f0); shading interp

hold on
contour(squeeze(Jbterm(:,:,2,t)),30, 'k');
hold off
colorbar;
set(gca, 'clim', [-1 1]*2);
% set(gca, 'clim', [-1 1].*1e-11)


%%
Jbterm = JBz;

jvec = reshape(Jbterm(:,:,2,20:100), nx*ny*81, 1);
% zvec = reshape(mfk(:,:,1,20:100).^2.*Hbl(:,:,1,40:120), nx*ny*81, 1);
% zvec = reshape((f0+ZETA(:,:,1,20:100)).*mfk(:,:,1,20:100), nx*ny*81, 1);
% zvec = reshape(ZETA(:,:,2,20:100).^2, nx*ny*81, 1);
zvec = reshape(f0*0.06*2*Hbl(:,:,1,20:100).*gradb(:,:,2,20:100).^2./f0, nx*ny*81, 1);
scatter((zvec), jvec);
title(num2str(corr( zvec, jvec)))
grid on
%%
[n bins] = histc( zvec ,10);

full(mean(sparse(1:length(zvec), bins, jvec)))

%%
t=10;
xp =70;
zlim = 1:40;
pcolor( 1:ny,Z,  squeeze(JBz(xp,:,:,t)).'); shading interp
cl = get(gca, 'clim');
hold on
contour(1:ny, Z(zlim), squeeze(repmat(Hbl(xp,:,1,t), [1 1 length(zlim) 1]).*Vz(xp,:,zlim,t).*by(xp,:,zlim,t)).', 50); 
plot(1:ny, -squeeze(Hbl(1,:,1,t)), 'k');
quiver(1:ny, Z, squeeze(V(xp,:,:,t)).', squeeze(W(xp,:,:,t)).', .25);
hold off
set(gca, 'clim', cl);
