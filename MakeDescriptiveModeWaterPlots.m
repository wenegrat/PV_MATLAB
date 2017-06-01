%%
% Load and Plot ARGO Climatology
filename = '~/Documents/RG_ArgoClim_Temperature_2016.nc';
TEMP = ncread(filename, 'ARGO_TEMPERATURE_MEAN');
filename = '~/Documents/RG_ArgoClim_Salinity_2016.nc';
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