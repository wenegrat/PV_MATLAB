%%%%%%%%%%%%%%%%% 

%Params
gap = [0.05 0.1]; margh=0.15; margw=0.1;
ti = 4*24/2; % GS_025_2F
ti = 191;
% ti = 200;
xi = 50; % GS_025_2F
xi = 36*2;
% xi = 130;
cl = [-10 10]; % GS_025_2F
cl = [-100 100];
zl = [-180 0];
xl = [35 65]; %GS_025_2F
xl = [60 72];
qon = true;
dc = 20; dcz = 5;
vf = 1e4; wf=vf;
mhs = 1e-4;
conts = [-1000 linspace(cl(1), cl(end), 19) 1000];
Bo = 9.81*2e-4*25./(1035*3994); % Adjust this depending on which model run comparing to.
eps = 0.1;
qs = 0.0;

norm = f0.*Bo./150;


% Form a slice line.
lc = [37 58];
rc = [46 67];
nl = 500;
xI = linspace(lc(1), rc(1), nl);
yI = linspace(lc(2), rc(2), nl);
distI = sqrt(xI.^2 + yI.^2);
distI = distI-distI(1);
distI = distI.*1000;
[XI, YI, ZI] = meshgrid(xI, yI, outputFull.Z);
theta = atan2(yI(end)-yI(1), xI(end)-xI(1));

JFi = NaN(nl, 50);
JDi = JFi;
Ti = JFi;
Vi = JFi;
Wi = JFi;
Qi = JFi;
Uprime = squeeze(U(:,:,:,ti)); - repmat(nanmean(U(:,:,:,ti)), [nx 1 1]);
Vprime = squeeze(V(:,:,:,ti)); - repmat(nanmean(V(:,:,:,ti)), [nx 1 1]);

for i=1:length(outputFull.Z)
JFi(:,i) = interp3(outputFull.X./1000, outputFull.Y./1000, outputFull.Z, permute(squeeze(outputFull.JFz(:,:,:,ti)),[2 1 3]), xI,yI,repmat(outputFull.Z(i), [1 nl]));
JDi(:,i) = interp3(outputFull.X./1000, outputFull.Y./1000, outputFull.Z, permute(squeeze(outputFull.JBz(:,:,:,ti)),[2 1 3]), xI,yI,repmat(outputFull.Z(i), [1 nl]));
 Ti(:,i) = interp3(outputFull.X./1000, outputFull.Y./1000, outputFull.Z, permute(squeeze(outputFull.T(:,:,:,ti)),[2 1 3]), xI,yI,repmat(outputFull.Z(i), [1 nl]));
 Vi(:,i) = interp3(outputFull.X./1000, outputFull.Y./1000, outputFull.Z, permute(squeeze(Vprime.*cos(pi/2-theta) + Uprime.*cos(theta)),[2 1 3]), xI,yI,repmat(outputFull.Z(i), [1 nl]));
 Wi(:,i) = interp3(outputFull.X./1000, outputFull.Y./1000, outputFull.Z, permute(squeeze(W(:,:,:,ti)),[2 1 3]), xI,yI,repmat(outputFull.Z(i), [1 nl]));
 Qi(:,i) = interp3(outputFull.X./1000, outputFull.Y./1000, outputFull.Z, permute(squeeze(outputFull.Q(:,:,:,ti)),[2 1 3]), xI,yI,repmat(outputFull.Z(i), [1 nl]));
end
%%
[nx, ny, nz, nt] = size(outputFull.Q);
% Plots

figure
%%%%% Plot 1 Surface D Flux
subtightplot(2, 3, 1, gap, margh, margw)

tprime = squeeze(outputFull.T(:,:,1,ti) - repmat(nanmean(outputFull.T(:,:,1,ti),1), [nx 1 1 1]));

[c, h]= contourf(outputFull.X./1000+eps, outputFull.Y./1000+eps, -squeeze(outputFull.JBz(:,:,2,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
hold on;
[c, h]= contourf(outputFull.X./1000, outputFull.Y./1000, -squeeze(outputFull.JBz(:,:,2,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');

yt = linspace(xl(1), xl(end), 100);
plot(xI, yI, 'k', 'LineWidth', 2);
%  contour(outputFull.X./1000, outputFull.Y./1000, squeeze(STRAIN(:,:,1,ti)).', [0 0], 'k', 'LineWidth', 1.);
hold off;
% set(gca, 'clim', [15 15.5]);
set(gca, 'clim', cl);

xlabel('x (km)'); 
ylabel('y (km)');
set(gca, 'FontSize', fs);
axis equal
title('$-\widehat{J_D}|_{z=0}$')
grid on

%%%%%%% PLOT 4 Slice D Flux
subtightplot(2,3, 4, gap, margh, margw);
[c, h]= contourf(distI, outputFull.Z+eps, -JDi.'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
hold on
[c, h]= contourf(distI, outputFull.Z, -JDi.'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
contour(distI, outputFull.Z, Ti.', 30, 'k', 'LineWidth', 1.25);
% contour(distI, outputFull.Z, Qi.', 30, 'k', 'LineWidth', 1.25);
set(gca, 'ylim', zl, 'xlim', [distI(1) distI(end)])
if (qon); quiver(distI(1:dc:end), outputFull.Z(1:dcz:end),vf.* squeeze(Vi( 1:dc:end, 1:dcz:end,1)).', wf.*squeeze(Wi(1:dc:end,1:dcz:end)).', qs,'k', 'LineWidth', 1., 'MaxHeadSize', mhs); end
% AnnotationQuiver(distI(1:dc:end), outputFull.Z(1:dcz:end), Vi(1:dc:end,1:dcz:end), Wi(1:dc:end,1:dcz:end), 100, 1, 1)
hold off
set(gca, 'clim', cl);
set(gca, 'XTick', (0:3:12).*1000, 'XTickLabel', strread(num2str(0:3:12),'%s'));
xlabel('Dist. (km)');
ylabel('z (m)');
set(gca, 'FontSize', fs);
grid on

%%%%% Plot 2 Surface F Flux
subtightplot(2, 3, 2, gap, margh, margw)

tprime = squeeze(outputFull.T(:,:,1,ti) - repmat(nanmean(outputFull.T(:,:,1,ti),1), [nx 1 1 1]));

[c, h]= contourf(outputFull.X./1000+eps, outputFull.Y./1000+eps, -squeeze(outputFull.JFz(:,:,2,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
hold on;
[c, h]= contourf(outputFull.X./1000, outputFull.Y./1000, -squeeze(outputFull.JFz(:,:,2,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');

yt = linspace(xl(1), xl(end), 100);
plot(xI, yI, 'k', 'LineWidth', 2);
%  contour(outputFull.X./1000, outputFull.Y./1000, squeeze(STRAIN(:,:,1,ti)).', [0 0], 'k', 'LineWidth', 1.);
hold off;
% set(gca, 'clim', [15 15.5]);
set(gca, 'clim', cl);

xlabel('x (km)'); 
ylabel('y (km)');
set(gca, 'FontSize', fs);
axis equal
title('$-\widehat{J_F}|_{z=0}$')
grid on

%%%%%%% PLOT 5 Slice F Flux
subtightplot(2,3, 5, gap, margh, margw);
[c, h]= contourf(distI, outputFull.Z+eps, -JFi.'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
hold on
[c, h]= contourf(distI, outputFull.Z, -JFi.'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
contour(distI, outputFull.Z, Ti.', 30, 'k', 'LineWidth', 1.25);
set(gca, 'ylim', zl, 'xlim', [distI(1) distI(end)])
if (qon); quiver(distI(1:dc:end), outputFull.Z(1:dcz:end),vf.* squeeze(Vi( 1:dc:end, 1:dcz:end,1)).', wf.*squeeze(Wi(1:dc:end,1:dcz:end)).', qs,'k', 'LineWidth', 1., 'MaxHeadSize', mhs); end
hold off
set(gca, 'clim', cl);
set(gca, 'XTick', (0:3:12).*1000, 'XTickLabel', strread(num2str(0:3:12),'%s'));
xlabel('Dist. (km)');
ylabel('z (m)');
set(gca, 'FontSize', fs);
grid on

%%%%% Plot 3 Surface Combined Flux
subtightplot(2, 3, 3, gap, margh, margw)

tprime = squeeze(outputFull.T(:,:,1,ti) - repmat(nanmean(outputFull.T(:,:,1,ti),1), [nx 1 1 1]));

[c, h]= contourf(outputFull.X./1000+eps, outputFull.Y./1000+eps, -squeeze(outputFull.JBz(:,:,2,ti)+outputFull.JFz(:,:,2,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
hold on;
[c, h]= contourf(outputFull.X./1000, outputFull.Y./1000, -squeeze(outputFull.JBz(:,:,2,ti)+outputFull.JFz(:,:,2,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');

yt = linspace(xl(1), xl(end), 100);
plot(xI, yI, 'k', 'LineWidth', 2);
%  contour(outputFull.X./1000, outputFull.Y./1000, squeeze(STRAIN(:,:,1,ti)).', [0 0], 'k', 'LineWidth', 1.);
hold off;
% set(gca, 'clim', [15 15.5]);
set(gca, 'clim', cl);

xlabel('x (km)'); 
ylabel('y (km)');
set(gca, 'FontSize', fs);
axis equal
title('$-(\widehat{J_F}+\widehat{J_D})|_{z=0}$')
grid on

%%%%%%% PLOT 6 Slice Combo Flux
subtightplot(2,3, 6, gap, margh, margw);
[c, h]= contourf(distI, outputFull.Z+eps, -(JFi+JDi).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
hold on
[c, h]= contourf(distI, outputFull.Z, -(JFi+JDi).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
contour(distI, outputFull.Z, Ti.', 30, 'k', 'LineWidth', 1.25);
set(gca, 'ylim', zl, 'xlim', [distI(1) distI(end)])
if (qon); quiver(distI(1:dc:end), outputFull.Z(1:dcz:end),vf.* squeeze(Vi( 1:dc:end, 1:dcz:end,1)).', wf.*squeeze(Wi(1:dc:end,1:dcz:end)).', qs,'k', 'LineWidth', 1., 'MaxHeadSize', mhs); end
hold off
set(gca, 'clim', cl);
set(gca, 'XTick', (0:3:12).*1000, 'XTickLabel', strread(num2str(0:3:12),'%s'));
cb= colorbar;
set(cb, 'Location', 'SouthOutside');
% set(get(cb, 'ylabel'), 'string', '$\frac{{-J}}{|fB_o/H_o|}$', 'Interpreter', 'Latex', 'Rotation', 0, 'FontSize', 20);
set(gca, 'clim', cl);
grid on

% colormap(cptcmap('cool-warm.cpt'));
colormap(cptcmap('balance.cpt'));

xlabel('Dist. (km)');
ylabel('z (m)');
set(gca, 'FontSize', fs);
set(gcf, 'Color', 'w', 'Position', [  665    44   857   924]);
% title('Slice at t = 3 days')
set(cb, 'Position', [ 0.3273    0.0441    0.3501    0.0260]);

