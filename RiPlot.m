
%% RiPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Plot of Ri
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Following Fox-Kemper Et al. 2008 definitions

by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);

%%
FrontalZoneMask = NaN(nx, ny, nz, nt);
FrontalZoneMask(:, [1:3 39:43 77:end], :, :) = 1; %Need to adjust this for a given domain.
deltaCrit = 0.01./(-2e-4*1035); % Convert to a delta T criteria.
tprime = THETA - repmat(THETA(:,:,1,:), [1, 1, nz, 1]);
mask = tprime<deltaCrit;
FrontalZoneMask(mask) = NaN;

bya = squeeze(nanmean(nanmean(nanmean(abs(by.*FrontalZoneMask)))));
bza = squeeze(nanmean(nanmean(nanmean(bz.*FrontalZoneMask))));


Ri = bza.*(1e-4).^2./bya.^2;
Rifull = bz.*(1e-4).^2./by.^2;
ShRed = (bya./1e-4).^2 - 4*bza;

%%
RiPlotFig = figure;
subplot(2,1,1)
plot(time, (Ri), 'LineWidth', 2)
set(gca, 'ylim', [-0 50], 'xlim', [0 25]) 
title('Ri_B')
set(gca, 'FontSize', fs)
grid on
subplot(2,1,2)
plot(time, bza, 'LineWidth',2 );
hold on
plot(time, (bya./(1e-4)).^2, 'LineWidth', 2); 
hold off
grid on
legend('N^2', 'Sh^2');
set(gca, 'FontSize', fs)


set(gcf, 'Color', 'w');