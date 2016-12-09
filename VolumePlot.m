zslice = [-1.51];
zslice = -5;
xslice = [0.25];
yslice = [0.25];
cl = [14 19];
cl = [15.25 18.1];
fs=14;
[x, y, z] = meshgrid(outputFull.X./1000, outputFull.Y./1000, outputFull.Z);
np = 4;
for i=1:np
    subplot(np/2,np/2,i)
tpos  = 1 + (i-1).*floor(5*86400./(7200));
% if i==np
%     tpos = 1 + floor(25*86400/7200);
% end
slice(x, y, z, permute(squeeze(outputFull.T(:,:,:,tpos)),[2,1,3]), xslice, yslice, zslice);
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
title(['Day : ', num2str(time(tpos) - time(1))], 'FontSize', fs);
% camzoom(1.4)
% camproj perspective
if i==1;
hold on
% plot3(squeeze(U(1,:,1,1)).*1e2 + 20, outputFull.Y./1000, 0*outputFull.Y./1000, 'k');
mask = NaN(150, 200);
mask(50,1:10:end) = 1;
% mask(30,2:20:end) =NaN;
quiver3(x(:,:,1), y(:,:,1), z(:,:,1)+1.5, repmat(squeeze(U(:,:,1,1)).*mask, [1 1]).', 0.*repmat(squeeze(U(:,:,1,1)), [1 1]).', 0.*repmat(squeeze(U(:,:,1,1)), [1 1]).', 30,...
    'lineWidth', 1.5, 'Color', 'k');
% quiver3(x, y, z, permute(U(:,:,:,tpos), [2, 1, 3]), permute(V(:,:,:,tpos), [2, 1, 3]), 0*permute(U(:,:,:,tpos), [2, 1, 3]));
hold off
end
set(gca, 'clim', cl)
box on
set(gca, 'BoxStyle', 'full');
end
% set(gca, 'zdir', 'reverse')
% view(-120, 25);
set(gcf, 'Color', 'w');

colormap(cptcmap('balance.cpt'))
