path = '/groups/thomas1/jacob13/GulfStream/NESEA/HIS/nesea_his.0062.nc';

W = GetVarROMS(path, path, {'omega', '(1)'}, {0, 0, [35 50], 0});
U = GetVarROMS(path, path, {'u', '(1)'}, {0, 0, [35 50], 0});
V = GetVarROMS(path, path, {'v', '(1)'}, {0, 0, [35 50], 0});

Eta = GetVarROMS(path, 0, {'zeta', '(1)'}, {0, 0, [1 1], 0});
theta_s = ncreadatt(path, '/', 'theta_s');
theta_b = ncreadatt(path, '/', 'theta_b');
hc = ncreadatt(path, '/', 'hc');
pm = ncread('/groups/thomas1/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc', 'pm'); %XXXXXXX
pn = ncread('/groups/thomas1/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc', 'pn'); %XXXXXXX

h = ncread('/groups/thomas1/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc', 'h'); %XXXXXXXXX

ti = 2;

z = compZ(path, 0, squeeze(Eta(:,:,ti)), theta_s, theta_b, hc, h);
Zx = DrvROMS(pm, z(:,:,35:end), 'x'); % XXX-Is this the right type of derivative to take?
Zy = DrvROMS(pn, z(:,:,35:end), 'y');

W = squeeze(W(:,:,:,ti)) + squeeze(U(:,:,:,ti)).*Zx + squeeze(V(:,:,:,ti)).*Zy;

[nx ny nz nt] = size(W);

W50 = NaN(nx, ny);
for i=1:nx
    disp(num2str(i))
    for j=1:ny
        W50(i,j) = interp1(squeeze(z(i,j,35:50)), squeeze(W(i,j,:)), -50);
    end
end

%%
wl = [-1 1].*600;
pcolor( (0:2001)./2,(0:1601)./2, squeeze(W50*86400).'); shading interp
hold on
% contour((0:2001)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
hold off
set(gca, 'clim',wl);
% title('$Temperature$');
set(gca, 'FontSize', 16);
axis equal
xlabel('km'); ylabel('km');
set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 
colormap( cptcmap('BlWhRe.cpt'))

title('Vertical velocity at 50 m depth', 'Interpreter', 'Latex', 'FontSize', 24);
cb = colorbar;
set(gcf, 'Color', 'w', 'Position', [           170         245        1112         719]);
set(cb, 'Position',[      0.90848201438849         0.157162726008345      0.0240719424460429         0.719054242002782]);
set(get(cb, 'ylabel'), 'string', 'm day$^{-1}$', 'Interpreter','latex', 'FontSize', 24)
% set(cb, 'Ticks', [0 100]);
set(gca, 'ylim', [0 600]);