%%%%%%%%%%%%%%%%% 

%Params
gap = [0.12 0.1]; margh=0.2; margw=0.1;
ti = 4*24/2; % GS_025_2F
ti = 30;
xi = 50; % GS_025_2F
xi = 130;
cl = [-10 10]; % GS_025_2F
cl = [-60 60];
zl = [-180 0];
xl = [35 65]; %GS_025_2F
xl = [10 90];
conts = [-100 linspace(cl(1), cl(end), 19) 100];
Bo = 9.81*2e-4*25./(1035*3994); % Adjust this depending on which model run comparing to.
eps = 0.1;

norm = f0.*Bo./150;

[nx, ny, nz, nt] = size(outputFull.Q);
% Plots

figure
% Plot 1 J_D vectors
subtightplot(3, 2, 1:2, gap, margh, margw)
plot(output.time, -(output.dJbsa+output.dJbea), 'LineWidth', 2, 'Color', 'k');
hold on
plot(output.time, -output.dJbsa, 'LineWidth', 2, 'LineStyle', '-.', 'Color', 'k');
plot(output.time, -output.dJbea,'LineWidth', 2, 'LineStyle', '--', 'Color', 'k');
yt = get(gca, 'YTick');
plot(output.time(ti).*ones(size(yt)), yt, 'k');
hold off
grid on
set(gca, 'xlim', [0 output.time(end)]);
l = legend('$-\int_A J^z_D$', '$-\int_A J^z_{D_{SURF}}$', '$-\int_A J^z_{D_{EDDY}}$', 'Location', 'SouthEast');
set(l, 'FontSize', fs);
xlabel('Day');
title('Diabatic PV Flux');
ylabel('$m^3s^{-4}$'); 
set(gca, 'FontSize', fs);
set(gca, 'ylim', [-3 0].*1e-3);
subtightplot(3, 2, [3 5], gap, margh, margw);
tprime = squeeze(outputFull.T(:,:,1,ti) - repmat(nanmean(outputFull.T(:,:,1,ti),1), [nx 1 1 1]));

[c, h]= contourf(outputFull.X./1000+eps, outputFull.Y./1000+eps, -squeeze(outputFull.JBz(:,:,2,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
hold on;
[c, h]= contourf(outputFull.X./1000, outputFull.Y./1000, -squeeze(outputFull.JBz(:,:,2,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');

yt = linspace(xl(1), xl(end), 100);
plot(outputFull.X(xi)./1000.*ones(size(yt)), yt, 'k', 'LineWidth', 2);

%  contour(outputFull.X./1000, outputFull.Y./1000, squeeze(outputFull.T(:,:,1,ti)).', 20, 'k', 'LineWidth', 1.25);
%  contour(outputFull.X./1000, outputFull.Y./1000, squeeze(ZETA(:,:,1,ti)).', 8, 'k', 'LineWidth', 1.);

hold off;
% set(gca, 'clim', [15 15.5]);
set(gca, 'clim', cl);

xlabel('x (km)'); 
ylabel('y (km)');
set(gca, 'FontSize', fs);

% PLOT 3
subtightplot(3,2, [4 6], gap, margh, margw);
vprime = V(xi,:,:,ti); - repmat(nanmean(V(xi,:,:,ti),2), [1 ny 1 1]);

% subplot(1,2,1)
[c, h]= contourf(outputFull.Y./1000+eps, outputFull.Z+eps, -squeeze(outputFull.JBz(xi,:,:,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');

% cl = get(gca, 'clim');
hold on
[c, h]= contourf(outputFull.Y./1000, outputFull.Z, -squeeze(outputFull.JBz(xi,:,:,ti)).'./norm, conts); shading interp
set(h, 'edgecolor', 'none');
 contour(outputFull.Y./1000, outputFull.Z, squeeze(outputFull.T(xi,:,:,ti)).', 30, 'k', 'LineWidth', 1.25);
dc = 4; dcz = 2;
set(gca, 'ylim', zl, 'xlim', xl)


quiver(outputFull.Y(1:dc:end)./1000, outputFull.Z(1:dcz:end), squeeze(vprime(1, 1:dc:end, 1:dcz:end,1)).', 20*squeeze(W(xi,1:dc:end,1:dcz:end,ti)).', 0.75,'k', 'LineWidth', 1.5)
% ncquiverref(Y(1:dc:end)./1000, Z(1:dcz:end), 100*squeeze(vprime(1, 1:dc:end, 1:dcz:end,1)).', 1000*squeeze(W(xi,1:dc:end,1:dcz:end,ti)).', '$cm/s V, x 10^{-3} cm/s W$',1,true,'k')
hold off
cb= colorbar;
set(cb, 'Location', 'SouthOutside');
set(get(cb, 'ylabel'), 'string', '$\frac{{-J_D}}{|fB_o/H_o|}$', 'Interpreter', 'Latex', 'Rotation', 0, 'FontSize', 20);

set(gca, 'clim', cl);

% colormap(cptcmap('cool-warm.cpt'));
colormap(cptcmap('balance.cpt'));

xlabel('y (km)');
ylabel('z (m)');
set(gca, 'FontSize', fs);
set(gcf, 'Color', 'w', 'Position', [   670   245   633   726]);
% title('Slice at t = 3 days')
set(cb, 'Position', [0.3250    0.0981    0.3501    0.0260]);

