%% dQdt Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dQdt Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dQdtFig = figure;
subplot(3,1,1)
[c h] = contourf(time, Y/1000, squeeze(THETA(xpos,:,1,:)), 50); 
set(h, 'edgecolor', 'none')
hold on
contour(time, Y/1000, squeeze(THETA(xpos,:,1,:)), isoT, 'k', 'LineWidth', 2);
hold off


ylabel('y (km)');
xlabel('Days');
grid on
set(gca, 'FontSize', fs);

subplot(3,1,2:3)
plot(time, smooth(Qt,1), 'LineWidth', 2)
hold on
plot(time, -dJFdt, 'LineWidth', 2);
plot(time, -dJBdt, 'LineWidth', 2);


plot(time,-(dJFdt + dJBdt), 'LineWidth', 2, 'LineStyle', '--');
plot(time, - dJBzdt)
% plot(qdira);
% plot(Qta, 'LineWidth', 2)

hold off
legend('\partialQ/\partial t', 'Fric', 'Dia', 'Sum',  'Location', 'SouthEast');
xlabel('Days');
ylabel('\partialQ/\partial t');
grid on
set(gca, 'FontSize', fs);
% cb = colorbar;
% set(cb, 'visible', 'off')
set(gcf, 'Color', 'w', 'Position', [675   445   830   530]);