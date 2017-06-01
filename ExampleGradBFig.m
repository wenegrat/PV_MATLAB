%%
bx = TtoB.*DPeriodic(output.Tsurf, 500, 'x');
by = TtoB.*DPeriodic(output.Tsurf, 500, 'y');
magb = sqrt(bx.^2 + by.^2);
magbn = (magb./f0.^2);
%%
figure

%  colormap(cptcmap('balance.cpt'))
 colormap(cptcmap('seminf-haxby.cpt'))

 cl = [0 15]
 conts = linspace(cl(1), cl(2), 20);

 tcl = [14.8 15.4];
 tconts = linspace(tcl(1), tcl(2), 20);
fs = 18;
gap = [.05 .05]; margh = .1; margw=.1;
for i=1:3

% if (i==nf)
%     tpos = length(time);
 %end
 tpos = 1 + (i-1).*floor(7.5*86400./(7200))+floor(5*86400./(7200));

 
subtightplot(1,3,i, gap, margh, margw)
% [c h] = contourf(outputFull.X/1000, outputFull.Y/1000, squeeze(outputFull.T(:,:,1,tpos)).', conts); 
[c h] = contourf(outputFull.X/1000, outputFull.Y/1000, squeeze(magbn(:,:,tpos)).', conts); 

% if ~(i==1)
%     set(gca, 'YTickLabel', []);
% else

% end
set(h, 'edgecolor', 'none')
hold on
contour(outputFull.X/1000, outputFull.Y/1000, squeeze(output.Tsurf(:,:,tpos)).',15.1:.05:15.3, 'k', 'LineWidth', 1);
hold off
% ylabel('y (km)');
xlabel('x (km)');
grid on
set(gca, 'FontSize', fs);
if i==1
    ylabel('y (km)');
end
if i==2
title('$\frac{|\nabla_h b|}{f^2}$', 'FontSize', 26);
end

axis equal
set(gca, 'clim', cl);
% cb = colorbar;
%     set(get(cb, 'ylabel'), 'String', '$\frac{\zeta}{f}$', 'Rotation', 0, 'FontSize', 28,'Interpreter', 'Latex');

end
cb = colorbar;
set(cb, 'Position', [ 0.91930401136895         0.281468531468531  0.0218724592192848         0.437194308105647]);

    set(gcf,'Color','w', 'Position', [675   404   789   570])