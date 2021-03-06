%% QBudgetPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QBudgetPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cl = [15.8 16.85];
fs = 16;
xpos = 16;

QBudgFig =figure;
subplot(2,2,1:2)
plot(time, Qa, 'LineWidth', 3)
hold on
plot(time, -JFa, 'LineWidth', 2); 
plot(time, -JBa, 'LineWidth', 2);
plot(time, -(JFa+JBa), 'LineWidth',3, 'LineStyle', '--');

hold off
legend('Q', '-J_F^z', '-J_B^z','Sum', 'Location', 'SouthWest');
xlabel('Days');
ylabel('\Delta Q');
grid on
title(titleString)
set(gca, 'FontSize', fs);

subplot(2,2,3)
contourf(X/1000, Y/1000, squeeze(THETA(:,:,1,tind)).'); shading interp
set(gca, 'clim', cl)
xlabel('x (km)'); ylabel('y (km)');
hold on
contour(X/1000, Y/1000, squeeze(THETA(:,:,1,tind)).', isoT, 'k', 'LineWidth', 2)

axis equal
grid on
set(gca, 'xlim', [0 X(end)/1000], 'ylim', [0 Y(end)/1000]);
set(gca, 'FontSize', fs);
yt = Y/1000;
plot(xpos*ones(size(yt)), yt,'--k', 'LineWidth', 2)
hold off
title(['Day: ', num2str(time(tind), 2)])

subplot(2,2,4)
[c, h]=contourf(Y/1000, Z, squeeze(THETA(xpos,:,:,tind)).', 20); 
% set(h, 'edgecolor', 'none');
set(gca, 'clim', cl)
xlabel('y (km)'); ylabel('z (m)');
hold on
contour(Y/1000, Z, squeeze(THETA(xpos,:,:,tind)).', isoT, 'k', 'LineWidth', 2)
% set(gca, 'ydir', 'reverse')
cb = colorbar;
set(get(cb, 'ylabel'), 'string', 'T (\circ C)');
set(gca, 'FontSize', fs);
grid on
set(gcf, 'Renderer', 'painters');
set(gcf, 'Color', 'w', 'Position', [675   370   825   605]);



