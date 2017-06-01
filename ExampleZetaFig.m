i=1;
figure
 colormap(cptcmap('balance.cpt'))
 
 cl = [-1 1]
 conts = linspace(cl(1), cl(2), 20);

 tcl = [14.8 15.4];
 tconts = linspace(tcl(1), tcl(2), 20);
tpos = 1 + (i-1).*floor(15*86400./(7200))+floor(5*86400./(7200));

% if (i==nf)
%     tpos = length(time);
 %end
 
 subplot(1,2,1)
% [c h] = contourf(outputFull.X/1000, outputFull.Y/1000, squeeze(outputFull.T(:,:,1,tpos)).', conts); 
[c h] = contourf(outputFull.X/1000, outputFull.Y/1000, squeeze(output.Tsurf(:,:,tpos)).', tconts); 

% if ~(i==1)
%     set(gca, 'YTickLabel', []);
% else
    ylabel('y (km)');
% end
set(h, 'edgecolor', 'none')
hold on
contour(outputFull.X/1000, outputFull.Y/1000, squeeze(output.Tsurf(:,:,tpos)).',15.1:.05:15.3, 'k', 'LineWidth', 1);
hold off
% ylabel('y (km)');
xlabel('x (km)');
grid on
set(gca, 'FontSize', fs);
title('$T\;[^{\circ} C]$', 'FontSize', fs);
axis equal
set(gca, 'clim', tcl);
cb = colorbar;
%     set(get(cb, 'ylabel'), 'String', '$T\;[^{\circ} C]$', 'Rotation', 0, 'FontSize', 28,'Interpreter', 'Latex');
set(gca, 'FontSize', 18)

subplot(1,2,2)
% [c h] = contourf(outputFull.X/1000, outputFull.Y/1000, squeeze(outputFull.T(:,:,1,tpos)).', conts); 
[c h] = contourf(outputFull.X/1000, outputFull.Y/1000, squeeze(ZETA(:,:,1,tpos)./f0).', conts); 

% if ~(i==1)
%     set(gca, 'YTickLabel', []);
% else
    ylabel('y (km)');
% end
set(h, 'edgecolor', 'none')
hold on
contour(outputFull.X/1000, outputFull.Y/1000, squeeze(output.Tsurf(:,:,tpos)).',15.1:.05:15.3, 'k', 'LineWidth', 1);
hold off
% ylabel('y (km)');
xlabel('x (km)');
grid on
set(gca, 'FontSize', fs);
title('$\zeta/f$', 'FontSize', fs);
axis equal
set(gca, 'clim', cl);
cb = colorbar;
%     set(get(cb, 'ylabel'), 'String', '$\frac{\zeta}{f}$', 'Rotation', 0, 'FontSize', 28,'Interpreter', 'Latex');
set(gca, 'FontSize', 18)
    set(gcf,'Color','w', 'Position', [675   404   789   570])