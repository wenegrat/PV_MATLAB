
%% RiPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Plot of Ri
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Following Fox-Kemper Et al. 2008 definitions

by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);

%%
[nx ny nz nt] = size(by);
bym = nanmedian(reshape(abs(by(:,:,2,:)), nx*ny, nt));
bymm = repmat(bym, [nx 1 ny nz]);
bymm = permute(bymm, [1 3 4 2]);
%%
FrontalZoneMask = NaN(nx, ny, nz, nt);
% FrontalZoneMask(:, [1:3 39:43 77:end], :, :) = 1; %Need to adjust this for a given domain.
% FrontalZoneMask(:,:,1:40,:) = 1;
% FrontalZoneMask = masknan;

FrontalZoneMask = (abs(by) > 1.1*abs(bymm));
deltaCrit = 0.01./(-2e-4*1035); % Convert to a delta T criteria.
tprime = THETA - repmat(THETA(:,:,1,:), [1, 1, nz, 1]);
mask = tprime<deltaCrit;
% FrontalZoneMask(mask) = NaN;
FrontalZoneMask = FrontalZoneMask.*masknan;
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