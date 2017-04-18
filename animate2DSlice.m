for i=1:nt
ti = i;
subplot(3,1,1)
contourf(1:200, outputFull.Z, squeeze(outputFull.Q(2,:,:,ti)).');
hold on
contour(1:200, outputFull.Z, squeeze(outputFull.T(2,:,:,ti)).', 'k');
contour(1:200, outputFull.Z, squeeze(outputFull.Q(2,:,:,ti)).',[0 0], 'r');

set(gca, 'clim', [-1 1].*1e-8);
plot(1:200, -squeeze(Hkpp(2,:,1,ti)), 'k', 'LineWidth', 2)
hold off
title(num2str(i))

subplot(3,1,2)
[c, h] = contourf(1:200, outputFull.Z, squeeze(W(2,:,:,ti)).');
set(h, 'edgecolor', 'none')
hold on
contour(1:200, outputFull.Z, squeeze(outputFull.T(2,:,:,ti)).', 30,'k');
set(gca, 'clim', [-1 1].*1e-4);
plot(1:200, -squeeze(Hkpp(2,:,1,ti)), 'k', 'LineWidth', 2)
hold off


subplot(3,1,3)
%%
[c, h] = contourf(fzone, outputFull.Z, squeeze(bprime(:,:,ti).*wprime(:,:,ti)).', 30);
set(h, 'edgecolor', 'none')
hold on
contour(1:200, outputFull.Z, squeeze(outputFull.T(2,:,:,ti)).', 30,'k');
contour(1:200, outputFull.Z, squeeze(Rib(2,:,:,ti)).',[1 1], 'r');

dz = 1;
dx = 3;
quiver(1:dx:200, outputFull.Z(1:dz:end), squeeze(V(2,1:dx:end,1:dz:end,ti)).', 500*squeeze(W(2,1:dx:end, 1:dz:end,ti)).',2,...
    'ShowArrowhead', 'off', 'LineWidth', 2, 'Color','r');

plot(1:200, -squeeze(Hkpp(2,:,1,ti)), 'k', 'LineWidth', 2)

hold off
% axis equal
set(gca, 'clim', [-1 1].*1e-6, 'xlim', [75 125], 'ylim', [-150 -0]);
pause(.01)
%%
end


%%
figure
subplot(2,1,1)

contourf(1:nt, Zl, squeeze(nanmean(outputFull.Q(2,fzone,:,:))),50); shading interp
colorbar;
hold on
contour(1:nt, Zl, squeeze(nanmean(outputFull.Q(2,fzone,:,:))),[ 0 0],'r'); shading interp
plot(1:nt, -squeeze(nanmean(Hkpp(2,fzone,1,1:end-1))),'k', 'LineWidth', 2)
hold off
set(gca, 'clim', [-1 1].*1e-9)
title('q');

subplot(2,1,2)
% [c, h] = contourf(1:nt, Zl, bpwp, 100);
% set(h, 'edgecolor', 'none')
% set(gca, 'clim', [-1 1].*1e-7)
% [c, h] = contourf(1:nt, Zl, squeeze(nanmean(outputFull.JBz(2,fzone,:,:)+outputFull.JFz(2,fzone,:,:))), 300);
% set(h, 'edgecolor', 'none')
% set(gca, 'clim', [-1 1].*5e-14)
% [c, h] = contourf(1:nt, Zl, log10(squeeze(nanmean(bz(2,fzone,:,:)))), 300);
% set(h, 'edgecolor', 'none')
% set(gca, 'clim', [-6 -4])
% [c, h] = contourf(1:nt, Zl, real(squeeze(nanmean(sqrt(f0.*outputFull.Q(2,fzone,:,:)./bz(2,fzone,:,:)))))./f0, 300);
% set(h, 'edgecolor', 'none')
% set(gca, 'clim', [-1 1].*1)
% [c, h] = contourf(1:nt, Zl, squeeze(nanmean(Vz(2,fzone,:,:)))./((6*f0).^2./f0), 300);
% set(h, 'edgecolor', 'none')
% set(gca, 'clim', [-1 1].*.5)
% [c, h] = contourf(1:nt, Zl, squeeze(bpwp), 300);
% set(h, 'edgecolor', 'none')
% set(gca, 'clim', [-1 1].*1e-7)
% [c, h] = contourf(1:nt, Zl, squeeze(gradient(bpwp.', Zl, 1e10).'), 300);
% set(h, 'edgecolor', 'none')
% set(gca, 'clim', [-1 1].*1e-8)
[c, h] = contourf(1:nt, Zl, squeeze(gradient(gradient(bpwp.', Zl, 1e10),Zl,1e10).'), 300);
set(h, 'edgecolor', 'none')
set(gca, 'clim', [-1 1].*1e-10)


hold on
contour(1:nt, Zl, squeeze(nanmean(outputFull.Q(2,fzone,:,:))),[ 0 0],'r'); shading interp

plot(1:nt, -squeeze(nanmean(Hkpp(2,fzone,1,1:end-1))),'k', 'LineWidth', 2)
hold off
colorbar
title('GSP')