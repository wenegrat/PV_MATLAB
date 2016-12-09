%% JOINT DIABATIC EDDY FLUX PLOTS
gap = [0.06 0.1]; margh = 0.15; margw=0.1;

figure
subtightplot(2, 5, [1 2 6 7], gap, margh, margw)
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
set(gca, 'ylim', zl, 'xlim', xl)
quiver(Y(1:dc:end)./1000, Z(1:dcz:end), squeeze(vprime(1, 1:dc:end, 1:dcz:end,1)).', 100*squeeze(W(xi,1:dc:end,1:dcz:end,ti)).', 0.5,'k', 'LineWidth', 1.5)
% ncquiverref(Y(1:dc:end)./1000, Z(1:dcz:end), 100*squeeze(vprime(1, 1:dc:end, 1:dcz:end,1)).', 1000*squeeze(W(xi,1:dc:end,1:dcz:end,ti)).', '$cm/s V, x 10^{-3} cm/s W$',1,true,'k')
hold off
cb= colorbar;
set(get(cb, 'ylabel'), 'string', '$\frac{{-J_D}}{f|B_o|/H}$', 'Interpreter', 'Latex', 'Rotation', 0, 'FontSize', 28);

set(gca, 'clim', cl);

colormap(cptcmap('cool-warm.cpt'));
xlabel('y (km)');
ylabel('z (m)');
set(gca, 'FontSize', fs);
set(gcf, 'Color', 'w');
title('Slice at t = 3 days')


% LINE PLOTS
xl = [0 30];



subtightplot(2,5,3:5, gap, margh, margw)
plot(time, output.dJba_t, 'LineWidth', 2);
hold on
plot(time, output.dJbsa, 'LineWidth', 2);
plot(time, output.dJbea,'LineWidth', 2);
% plot(time, output.dJb, 'LineWidth', 2);
hold off
l=legend('$-\int J_{D_{Total}}\mathrm{dA}$', '$-\int J_{D_{SURF}}\mathrm{dA}$', '$-\int J_{D_{EDDY}}\mathrm{dA}$', 'Location', 'NorthWest');
set(l, 'Interpreter', 'Latex', 'FontSize', 20);
grid on

set(gca, 'xlim', xl);
set(gca, 'XTickLabel', []);
set(gca, 'FontSize', fs)
ylabel('$m^{3}$ $s^{-4}$', 'FontSize', 19);
% title('Integrated over surface')

subtightplot(2,5,8:10, gap, margh, margw)
KEd = -(KE-KE(1));
APEd = -(APE.'-APE(1));
% plot(time, gradient(KEd, ts),'LineWidth', 2);
% hold on
% plot(time, gradient(APEd, ts),'LineWidth', 2);
plot(time, KEd.*dx.*dy,'LineWidth', 2);
hold on
plot(time, APEd.*dx.*dy,'LineWidth', 2);
% plot(time, -(KE-KE(1))-(APE.'-APE(1)));
% plot(time, PE,'LineWidth', 2);
% plot(time, cumtrapz(SC)*ts);
% plot(time, APEm-APEm(1) + (KEm-KEm(1)))
hold off
grid on
% l = legend('$\frac{\partial KE}{\partial t}$', '$\frac{\partial APE}{\partial t}$');
l = legend('$\Delta \int KE \,\mathrm{dV}$', '$\Delta \int APE\,\mathrm{dV}$', 'Location', 'NorthWest');

set(l, 'Interpreter', 'Latex', 'FontSize', 19);
set(gca, 'xlim', xl);
set(gca, 'FontSize', fs)
xlabel('Days', 'FontSize', 20);
ylabel('$kg\,m^2$ $s^{-2}$', 'FontSize', 19);
set(gcf, 'Color', 'w', 'Position', [  289         310        1495         607])