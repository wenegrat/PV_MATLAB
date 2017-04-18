%% QBudgetPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QBudgetPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
op = outputFull;
[nx, ny, nz, nt] = size(THETA);
timelimit = [0 30];
cl = [13.75 16.75];
cl = [13.5 17];
fs = 16;
xpos = 16;
zslice = [-1.51];
zslice = -5;
xslice = [0.25];
yslice = [0.25];

dh = diff([Zl; -300]);
gridvol = permute(repmat( dx.*dy.*abs(dh), [1 nx ny nt]), [2 3 1 4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL VOlUME ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isoT = [13.75 15.75]; %Needed for plot below
% isoT = [13.5 15.5];
mask = double((THETA(:,:,:,:)>isoT(1)) & (THETA(:,:,:,:)<isoT(2)));
vol = squeeze(sum(sum(sum(mask.*gridvol))));

Q = op.Q;
JFz = op.JFz; JBz = op.JBz;
ts = op.ts;
dx = op.dx; dy = op.dy;

%INTEGRATE Q
IntegrateQTerms;
%AreaIntegrateJTerms
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JFb, ~    ] = areaIntegrateJVecs(squeeze(JFz(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JFa = JFs-JFb;
[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JBb, ~    ] = areaIntegrateJVecs(squeeze(JBz(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JBa = JBs-JBb;

% Surface Only Scaling
JS = -1e-4.*2e-4.*9.81./(3994*1035).*Qo(:,:,:,1:end-1)./H(:,:,:,1:end-1);
JS = 1.2.*JS; %Efficiency Factor
[JFS, ~    ] = areaIntegrateJVecs(squeeze(JS(:,:,end,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);

[x, y, z] = meshgrid(op.X./1000,op.Y./1000, op.Z);

QBudgFig =figure;
% colormap(cptcmap('balance.cpt'))

gap = [0.05 0.075]; margh=.1; margw=.15;
%% VOLUME PLOTS
for i=1:3
    subtightplot(5,3,[i i+3], gap, margh, margw)
tpos  = 1 + (i-1).*floor(5*86400./(7200));
tpos = 1 + (i-1).*floor(5*86400./(7200))+floor(5*86400./(7200));

% if i==np
%     tpos = 1 + floor(25*86400/7200);
% end
slice(x, y, z, permute(squeeze(THETA(:,:,:,tpos)),[2,1,3]), xslice, yslice, zslice);
if (i==1)
%     cl = get(gca, 'clim');
end
set(gca, 'clim', cl);
shading interp
set(gca, 'xlim', [0 80], 'ylim', [0 100], 'zlim', [-300 0]);
daspect([1,1,8])
axis tight
view(-16, 25);
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (m)');
title(['Day : ', num2str(op.time(tpos) - op.time(1))], 'FontSize', fs);
% camzoom(1.4)
% camproj perspective
hold on
h = contourslice(x, y, z, permute(squeeze(THETA(:,:,:,tpos)),[2,1,3]),xslice, yslice, zslice, isoT);
set(h, 'EdgeColor', [1 1 1].*0.35);
set(h, 'LineWidth', 2);
% quiver3(x, y, z, permute(U(:,:,:,tpos), [2, 1, 3]), permute(V(:,:,:,tpos), [2, 1, 3]), 0*permute(U(:,:,:,tpos), [2, 1, 3]));
hold off
set(gca, 'clim', cl)
box on
set(gca, 'BoxStyle', 'full');
end

cb = colorbar;
set(get(cb, 'ylabel'), 'String', 'T $[^{\circ}C]$', 'FontSize', 18, 'Interpreter', 'Latex', 'Rotation', 0)
set(cb, 'Ticks', [13.75 15.25 16.75]);
set(cb, 'ticklabels', {17 18.5 20});
set(cb, 'FontSize', 14);
%% Line PLOTS
subtightplot(5,3,7:9, gap, margh, margw)
plot(outputFull.time, -squeeze(Qo(1,1,1,1:end-1)), 'LineWidth', 2);
hold on
% plot(outputFull.time, -squeeze(Qsw(1,1,1,1:end-1)), 'LineWidth', 2);
% plot(outputFull.time, -(squeeze(Qo(1,1,1,1:end-1))-squeeze(Qsw(1,1,1,1:end-1))), 'LineWidth', 2);
hold off
grid on
set(gca, 'ylim', [0 125]);
% set(gca, 'ylim', [-250 250]);
ylabel('$Q_o$ $(W m^{-2})$');
set(gca, 'FontSize', fs);
set(gca, 'xlim', timelimit);
set(gca, 'XTickLabel', []);

subtightplot(5,3,10:15, gap, margh, margw)
plot(outputFull.time, Qa, 'LineWidth', 3)
hold on
plot(outputFull.time, -JFa, 'LineWidth', 2); 
plot(outputFull.time, -JBa, 'LineWidth', 2);
% plot(outputFull.time, -(JFa+JBa), 'LineWidth',3, 'LineStyle', '--');
plot(outputFull.time, -JFS, 'LineWidth', 1.5, 'LineStyle','--');

set(gca, 'xlim', timelimit);

hold off
lgd = legend('$\Delta q$', '$-\int_{t} J_F^z$', '$-\int_t J_D^z$','$-\int_t J_{D}^{{SURF}}$', 'Location', 'NorthEastOutside');
xlabel('Days');
ylabel('$\Delta q$ $(s^{-3})$');
grid on
set(gca, 'FontSize', fs);
hold on
hold off
set(gca, 'ylim', [-1 1].*2.5e-9)



set(gcf, 'Color','w', 'Position', [  324         270        1090         679]);

set(lgd, 'Position', [0.87   0.1954    0.0743    0.1848]);


% cmin = cptcmap('temperature(1).cpt');
% xi = 1:1:length(cmin);
% xo = linspace(1, length(cmin), 100);
% for i=1:3
%     cout(:,i) = interp1(xi, cmin(:,i), xo);
% end
% colormap(cout);
colormap(cptcmap('balance.cpt'));
set(cb, 'Position', [ 0.8779    0.6156    0.0148    0.2900]);
