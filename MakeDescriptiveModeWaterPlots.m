%%
% Load and Plot ARGO Climatology
filename = '/data/thomas/jacob13/RG_ArgoClim_Temperature_2016.nc';
TEMP = ncread(filename, 'ARGO_TEMPERATURE_MEAN');
filename = '/data/thomas/jacob13/RG_ArgoClim_Salinity_2016.nc';
SAL = ncread(filename, 'ARGO_SALINITY_MEAN');
TIME = ncread(filename, 'TIME');

[nlon, nlat, np] = size(TEMP);
% TA = NaN(nlat, nlon, 12);
% for i=1:12
%     inds = 1:12:(nt-12) + (i-1);
%    TA(:,:,i) = nanmean(TEMP(:,:,inds), 3);
% end
PRESS = ncread(filename, 'PRESSURE');
LAT = ncread(filename, 'LATITUDE');
LON = ncread(filename, 'LONGITUDE');

%%
SA = NaN(nlon, nlat, np);
CT = SA;
RHO = SA;
for i=1:nlat
    for j=1:nlon
SA(j,i,:) = gsw_SA_from_SP(squeeze(SAL(j,i,:)), PRESS,mod( LON(j), 360), LAT(i));
CT(j,i,:) = gsw_CT_from_t(squeeze(SA(j,i,:)), squeeze(TEMP(j,i,:)), PRESS);
RHO(j,i,:) = gsw_rho(squeeze(SA(j,i,:)), squeeze(CT(j,i,:)), 0);
    end
end
%%
dz = diff(PRESS);
dsdz = NaN(nlon, nlat, np);
for i=1:nlat
    for j=1:nlon
        
        dsdz(j,i,1:end-1) = -(squeeze(RHO(j,i,1:end-1)-RHO(j,i,2:end)))./dz;
    end
end
Q = 1./RHO.*dsdz.*permute(repmat( 2*2*pi./86400*sind(LAT), [1,nlon, np]), [2 1 3]);

%%
MT = NaN(nlon, nlat, 1);
for i=1:nlat;
    for j=1:nlon;
       indu = find(squeeze(TEMP(j,i,:))<=19,1);
       indb = find(squeeze(TEMP(j,i,:))<=17,1);
       if (~isempty(indu)&&~isempty(indb))
           MT(j,i,:) = PRESS(indb)-PRESS(indu);
       end
    end
end
    %%
filt = ones(3);
filt(2,2) = 2;
filt = filt./sum(sum(filt));
MTS = filter2(filt, MT);
%%

mask = ones(size(Q));
% mask(:,84:108,:)=1;
% mask = 1;
ind = find(LON==180.5,1);
varp =Q;
d = 25;
m_proj('mercator', 'longitudes', [-90 0], 'latitudes', [15 45]);
m_coast('patch',[.7 .7 .7]);
m_grid('box', 'fancy', 'fontsize', 26);
% m_contourf([LON; LON-360], LAT, [squeeze(MT(:,:,1)).' squeeze(MT(:,:,1)).'], 10);

m_pcolor([LON; LON-360], LAT, [squeeze(MT(:,:,1)).' squeeze(MT(:,:,1)).']);
% m_pcolor([LON; LON-360], LAT, [squeeze(Q(:,:,25)).' squeeze(Q(:,:,25)).']);

cb = colorbar;
set(get(cb, 'ylabel'), 'string', 'EDW Thickness (m)');
hold on
pvconts = linspace(0, 1e-10, 10);
c = m_contour([LON; LON-360], LAT, [squeeze(mask(:,:,25).*Q(:,:,25)).' squeeze(mask(:,:,25,:).*Q(:,:,25)).'], pvconts, 'LineColor', 'k', 'LineStyle', '--');
% set(c, 'linecolor', 'k')
% clabel(c);
hold off
shading interp
m_coast('patch',[.7 .7 .7]);
m_grid('box', 'fancy');
set(gca, 'FontSize', 16);
title('ARGO Climatology 2004-2015');
set(gcf, 'Color', 'w');

%% FANCY MAP
%% LOAD SHORELINE
shorelines = gshhs('~/GSHH/gshhs_c.b');%_h for high, _c for crude, _i for intermediate, from https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
levels = [shorelines.Level];
land = (levels==1);
%%
% xl = 600:1602;
% yl = 1:550;
% lats = lat(xl,yl); lons = lon(xl, yl);
% lonf = [lons(1,:) lons(:,end).' fliplr(lons(end,:)) flipud(lons(:,1)).'];
% latf = [lats(1,:) lats(:,end).' fliplr(lats(end,:)) flipud(lats(:,1)).'];
mask = ones(length(LON), length(LAT));
mask(repmat(LAT.',[length(LON) 1]) < 1) = NaN;

gshow = true;
cl = 5e-12;
pvconts = linspace(1e-12, 2e-10, 10);
pvconts2 = [1e-10 1e-10]*1;
cl = 1e-11; % For no hkpp bl mask
mc = .9; %coast color
% gap = [.05 .01], margh = .1; margw=.1;
figure
% subtightplot(2,2,1, gap, margh, margw);
axesm('mercator','MapLatLimit', [0 50], 'MapLonLimit', [-90 0], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 15, 'MLabelLocation', 15,'MLabelParallel', 'South',...
    'PLineLocation', 15, 'MLineLocation', 15, 'Frame', 'on', 'FontSize', 16);
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(LAT(:,:)),squeeze(LON(:,:)),(squeeze(MTS(:,:,1))).'); shading interp
contourm(LAT, LON, squeeze(mask.*squeeze(Q(:,:,25))).', pvconts, '--k')
% contourm(LAT, LON, squeeze(squeeze(TEMP(:,:,1))).', 16:19, '-k', 'LineWidth',2)

hold on
% hf = fillm(latf, lonf, 'k');
% set(hf, 'FaceColor', 'none');
% set(hf, 'LineWidth', 1, 'LineStyle', '--');
hold off
if gshow; geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc); end
set(gca, 'clim', [0 400]);
% t= textm(45, 360-82.5,0, '$T$', 'FontSize', 20, 'Color', 'k', ...
%     'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');
title({'\makebox[4in][c]{North Atlantic Subtropical Mode Water Thickness}',...
    '\makebox[4in][c]{[$17^{\circ}-19^{\circ}$, ARGO 2004-2015]}'},...
    'Interpreter','latex', 'FontSize', 18)
% title({'North Atlantic Subtropical Mode Water Thickness','[$17^{\circ}-19^{\circ}$, ARGO 2004-2015]'});
colormap(cptcmap('seminf-haxby.cpt'));
cb = colorbar;
set(get(cb, 'ylabel'), 'String', '$m$', 'Interpreter', 'Latex','Rotation', 0,'FontSize', 22);
set(cb, 'FontSize', 18);
% set(gca, 'FontSize', 20);
set(gcf, 'Color', 'w', 'Position', [675   369   852   605]);