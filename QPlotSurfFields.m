%% Q Plot 2 Panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q Plot 2 Panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

QPlotFig = figure;
fs = 12;
nf = 3;
cl = [15 18];
cl = [16.35 17.5];
cl = [-2 2];
conts = linspace(cl(1), cl(2), 20);
 colormap(cptcmap('balance.cpt'))

% colormap(cptcmap('balance.cpt'))
for i=1:nf
subplot(2,nf,i)
tpos = 1 + (i-1).*floor(7.5*86400./(7200));
tpos = 1 + (i-1).*floor(15*86400./(7200))+floor(5*86400./(7200));

% if (i==nf)
%     tpos = length(time);
 %end

% [c h] = contourf(outputFull.X/1000, outputFull.Y/1000, squeeze(outputFull.T(:,:,1,tpos)).', conts); 
[c h] = contourf(outputFull.X/1000, outputFull.Y/1000, squeeze(ZETA(:,:,1,tpos)./f0).', conts); 

% if ~(i==1)
%     set(gca, 'YTickLabel', []);
% else
    ylabel('y (km)');
% end
set(h, 'edgecolor', 'none')
hold on
% contour(time, Y/1000, squeeze(THETA(xpos,:,1,:)), isoT, 'k', 'LineWidth', 2);
hold off
% ylabel('y (km)');
xlabel('x (km)');
grid on
set(gca, 'FontSize', fs);
title(['Day: ', num2str(time(tpos) - time(1))], 'FontSize', fs);
axis equal
set(gca, 'clim', cl);
if i==3
    cb = colorbar;
    set(get(cb, 'ylabel'), 'String', '$\frac{\zeta}{f}$', 'Rotation', 0, 'FontSize', 24,'Interpreter', 'Latex');
end
end

%%
subplot(2,nf,(nf+1):(2*nf))
plot(output.time, smooth(output.Qa,1), 'LineWidth', 2)
hold on
plot(output.time, -output.Jfa   , 'LineWidth', 2);
plot(output.time, -output.Jba, 'LineWidth', 2);


plot(output.time,-(output.Jfa+ output.Jba), 'LineWidth', 2, 'LineStyle', '--');
% plot(time, - dJBzdt)
% plot(qdira);
% plot(Qta, 'LineWidth', 2)

hold off
% lgd = legend('\Delta q ', '-\int_{t} J_F^z', '-\int_t J_D^z','-\int_t J_{D_{SURF}}', 'Location', 'NorthWest');
xlabel('Days');
ylabel('$\Delta q$ $(s^{-3})$');
grid on
set(gca, 'FontSize', fs);
% cb = colorbar;
% set(cb, 'visible', 'off')
set(gcf, 'Color', 'w', 'Position', [ 670   440   660   532]);
set(gca, 'xlim', [0 time(end)])
set(cb , 'Position', [.92 .588 0.02 0.335])






