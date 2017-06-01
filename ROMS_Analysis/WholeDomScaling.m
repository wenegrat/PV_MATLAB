%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SCALING MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pardir = '/data/thomas/jacob13/GULFZ/'; % 1.5 km domain
basepath = [pardir 'HIS/'];
path1 = [pardir 'gulfz_grd.nc'];
forcepath = [pardir 'gulfz_frc.nc'];
ntp = 5;
offset = 0;
nt = 345; % Only looking at 1 year (should add offset?)

files = dir([basepath,'*.nc']);

tf = ncread([basepath, files(1).name], 'ocean_time'); % Get the first timestamp

zl = 47:50; %Only need to consider the surface domain XX- What is the minimum?
nz = length(zl);
xl = 1:1602; % Full Domain
yl = 1:922;
slice =  {[xl(1) xl(end)], [yl(1) yl(end)], [zl(1) zl(end)], 0};

% Load constants/Initial Params
path = [basepath, files(1).name];
time = ncread(path, 'ocean_time');
% nt = length(time);
ts = time(2)-time(1);
modeltime = tf(1):ts:((nt-1)*ts+tf(1));

% Load grid basics
path = path1;
f = ncread(path, 'f');
f = f(xl, yl);
rho0 = 1027.4;
g = 9.81;
pm = ncread(path, 'pm');
pm = pm(xl, yl);
pn = ncread(path, 'pn');
pn = pn(xl,yl);

h = ncread(path, 'h');

latst = ncread(path1, 'lat_rho'); %XXXXXXX
lat = latst(xl, yl);
lonst = ncread(path1, 'lon_rho'); %XXXXXXXXXXX
lon = lonst(xl, yl);
% FULL DOMAIN
% path = [basepath, files(50).name];
% 
% T = GetVarROMS(path, 0, {'temp', '(1)'}, {0, 0,0, 0});
% 
% subplot(2,1,1)
% pcolor(squeeze(T(:,:,24,1)).'); shading interp
% hold on
% rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
% hold off
% colorbar;

%% Calculate SCALING terms
[nx ny] = size(pm);

% Preallocate
Bx = NaN(nx, ny, nz, nt);
By = Bx;
Tf = Bx;
Sf = Tf;
Rho = Bx;
hkpp = NaN(nx, ny, nt);


name = files(1).name;
path = [basepath, name];
% Sigma coordinate parameters    
theta_s = ncreadatt(path, '/', 'theta_s');
theta_b = ncreadatt(path, '/', 'theta_b');
hc = ncreadatt(path, '/', 'hc');

h = h(xl, yl);

for i=1:nt;
    % 1) Housekeeping
    disp([num2str(i), '/', num2str((length(files)-ntp)*ntp)]);
    fileind = ceil((i+offset)/ntp);
    name = files(fileind).name;
    path = [basepath, name];
%     sliceind = ntp - mod(i,ntp);
    sliceind = mod(i,ntp);
    if sliceind==0; sliceind=ntp;end
    sliceT = {slice{1}, slice{2}, slice{3},[sliceind sliceind]};
    
    % 2) Calculate z coordinates at this timestep.
    Eta = GetVarROMS(path, 0, {'zeta', '(1)'}, sliceT);

    hkpp(:,:,i) = GetVarROMS(path, 0, {'hbls', '(1)'}, sliceT);
    
    z = compZ(path, 0, Eta, theta_s, theta_b, hc, h);
    z = z(:,:,zl);
    zwt = compZ(path, 1, Eta,  theta_s, theta_b, hc, h);
    [~, ~, nz] = size(z);
    zm = squeeze(nanmean(nanmean(z)));
    
    % 5) Load Buoyancy Terms
    T = GetVarROMS(path, 0, {'temp', '(1)'}, sliceT);
     Tf(:,:,:,i) = T;
    S = GetVarROMS(path, 0, {'salt', '(1)'}, sliceT);
    Sf(:,:,:,i) = S;

    rho = rho_eos(T, S, 0); % CROCO function (checked for consistency 1/12/17)
    Rho(:,:,:,i) = rho;
    B = -g*rho./rho0; % In-Situ B

    %Calculate Gradients
    Bx(:,:,:,i) = DrvS(pm, z, B, 'x');
    By(:,:,:,i) = DrvS(pn, z, B, 'y');
    
    
end


%% SURFACE FORCING TERMS
Qor = ncread(forcepath,'shflux'); % Monthly Values W/m^2
SWR = ncread(forcepath, 'swrad');
EPr = ncread(forcepath,'swflux')./(100.*86400); %Freshwater flux m/s
tau_u = ncread(forcepath,'sustr'); %N/m^2;
tau_u = Int_varROMS(tau_u, [2 1], [1 1]);
tau_v = ncread(forcepath, 'svstr'); %N/m^2
tau_v = Int_varROMS(tau_v, [3 1], [1 1]);

Qor = Qor(xl, yl,:);
SWR = SWR(xl, yl,:);
EPr = EPr(xl, yl,:);
tau_u = tau_u(xl, yl,:);
tau_v = tau_v(xl, yl,:);
SST = ncread(forcepath, 'SST');
SST = SST(xl, yl,:);
Qor = repmat(Qor, [1 1 2]);
EPr = repmat(EPr, [1 1 2]);
tau_u = repmat(tau_u, [1 1 2]);
tau_v = repmat(tau_v, [1 1 2]);
SST = repmat(SST, [1 1 2]);
SWR = repmat(SWR, [1 1 2]);
ds = 86400;
climtime = ((360+15)*ds):30*ds:(360*2*ds+360*ds);

[nx, ny, ~] = size(Qor);

Qo = NaN(nx, ny, nt);
EP = Qo;
SW = Qo;
tx = Qo;
sst = Qo;
ty = Qo;

[X, Y, T] = ndgrid(xl, yl, climtime);
[Xm, Ym, mt] = ndgrid(xl, yl, modeltime);
QG = griddedInterpolant(X, Y, T, Qor, 'cubic', 'none');
EPG = griddedInterpolant(X, Y, T, EPr, 'cubic', 'none');
TXG = griddedInterpolant(X, Y, T, tau_u, 'cubic', 'none');
TYG = griddedInterpolant(X, Y, T, tau_v, 'cubic', 'none');
SSTG = griddedInterpolant(X, Y, T, SST,'cubic', 'none');
SWG = griddedInterpolant(X, Y, T, SWR, 'cubic', 'none');

Qo = QG(Xm, Ym, mt);
EP = EPG(Xm, Ym, mt);
tx = TXG(Xm, Ym, mt);
ty = TYG(Xm, Ym, mt);
sst = SSTG(Xm, Ym, mt);
SW = SWG(Xm, Ym, mt);
% Qor = ncread(forcepath,'shflux'); % Monthly Values W/m^2
% EPr = ncread(forcepath,'swflux')./(100.*86400); %Freshwater flux m/s
% tau_um = ncread(forcepath,'sustr'); %N/m^2;
% tau_vm = ncread(forcepath, 'svstr'); %N/m^2
% Qor = Qor(xl, yl,:);
% EPr = EPr(xl, yl,:);
% tau_u(xl(2:end),yl, :) = tau_um;%(xl, yl,:);
% tau_v(xl, yl(2:end),:) = tau_vm;%(xl, yl(2:end),:);
% SST = ncread(forcepath, 'SST');
% SST = SST(xl, yl,:);
% Qor = repmat(Qor, [1 1 2]);
% EPr = repmat(EPr, [1 1 2]);
% tau_u = repmat(tau_u, [1 1 2]);
% tau_v = repmat(tau_v, [1 1 2]);
% SST = repmat(SST, [1 1 2]);
% climtime = datenum(2012,1, 15:30:(360*2));
% 
% % disp('here')
% Qo = NaN(nx, ny, nt);
% EP = Qo;
% tx = Qo;
% sst = Qo;
% ty = Qo;
% for x = 1:nx;
%     for y = 1:ny
% Qo(x, y,:) = interp1(climtime, squeeze(Qor(x,y,:)), modeltime, 'pchip');
% EP(x,y,:) = interp1(climtime, squeeze(EPr(x,y,:)), modeltime, 'pchip');
% tx(x,y,:) = interp1(climtime, squeeze(tau_u(x,y,:)), modeltime, 'pchip');
% ty(x,y,:) = interp1(climtime, squeeze(tau_v(x,y,:)), modeltime, 'pchip');
% sst(x,y,:) = interp1(climtime, squeeze(SST(x,y,:)), modeltime, 'pchip');
%     end
% end

Tmean = squeeze(nanmean(nanmean(nanmean(Tf(:,:,end,:), 4))));
Smean = squeeze(nanmean(nanmean(nanmean(Sf(:,:,end,:), 4))));

%% GEN THEORY SCALINGS
del = 1;
%Spatially and temporally variable alpha, beta
rhot = rho_eos(squeeze(Tf(:,:,end,:)), squeeze(Sf(:,:,end,:)), 0);
alpha = 1./rhot.*(rho_eos(squeeze(Tf(:,:,end,:))+del, squeeze(Sf(:,:,end,:)), 0)-rho_eos(squeeze(Tf(:,:,end,:))-del, squeeze(Sf(:,:,end,:)), 0))./(2*del);
beta = 1./rhot.*(rho_eos(squeeze(Tf(:,:,end,:)), squeeze(Sf(:,:,end,:))+del, 0)-rho_eos(squeeze(Tf(:,:,end,:)), squeeze(Sf(:,:,end,:))-del, 0))./(2*del);

Cp = 3994; % XXX - Get ROMS Value
ssta =  squeeze(Tf(:,:,end,:))-sst; % SST anomaly for surf flux
Bo = g*alpha.*(Qo-30*ssta)./(rho0*Cp) +  g*beta.*EP.*squeeze(Sf(:,:,end,:));

H = hkpp; H(:,:,1) = H(:,:,2);

M2 = squeeze( nanmean(Bx(:,:,:,:),3).^2 + nanmean(By(:,:,:,:),3).^2);

tmag = abs(tx+1i*ty);

H(H<6) = 6;

JFT = -0.2*M2.*H;
% H(H<24) = 24;
Hm =H;
JBTs = 1.2.*repmat(f, [1 1 nt]).*Bo./Hm;
% JBTs(JBTs<0) = 0;
JBTe = 0.15.*M2.*H;
JBT = JBTs + JBTe;

% disp('here')

R = 0.05.*(H.^2.*M2)./(Bo.*repmat(f, [1 1 nt]));

JFW = -tx./(rho0*H).*squeeze(By(:,:,end-2,:)) + ty./(rho0*H).*squeeze(Bx(:,:,end-2,:));

%% DECIMATED SCALINGS
fil = fspecial('average', [7 7]);
Hd = NaN(nx, ny, nt);
M2F = Hd;
ByF = Hd;
BxF = Hd;
Tx = Hd;
Ty = Hd;
BoF = Hd;
RhoF = Hd;
Tfilt = Hd;
for i=1:nt
    i
Hd(:,:,i) = filter2(fil, squeeze(H(:,:,i)));
BxF(:,:,i) = filter2(fil, squeeze(nanmean(Bx(:,:,:,i),3)));
ByF(:,:,i) = filter2(fil, squeeze(nanmean(By(:,:,:,i),3)));
M2F(:,:,i) = squeeze(BxF(:,:,i).^2 + ByF(:,:,i).^2);
Ty(:,:,i) = filter2(fil, squeeze(ty(:,:,i)));
Tx(:,:,i) = filter2(fil, squeeze(tx(:,:,i)));
BoF(:,:,i) = filter2(fil, squeeze(Bo(:,:,i)));
RhoF(:,:,i) = filter2(fil, squeeze(Rho(:,:,end-1,i)));
Tfilt(:,:,i) = filter2(fil, squeeze(Tf(:,:,end,i)));
end

JBTs_FILT = 1.2.*repmat(f, [1 1 nt]).*BoF./Hd;
JBTe_FILT = 0.15.*M2F.*Hd;
JFT_FILT = -0.2.*M2F.*Hd;
JFW_FILT = -Tx./(rho0.*Hd).*ByF + Ty./(rho0.*Hd).*BxF;



%%
iso = 1026.1;
delt = 0.2; 
mask = squeeze((Rho(:,:,end-1,:) > iso-delt) & (Rho(:,:,end-1,:) < iso+delt));
% mask = mask &hkpp>100; % Turn on/off for masking deep boundary layers as in Maze et al. 2013
DIA = sum(JBT.*mask, 3);
DIAs = sum(JBTs.*mask,3);
DIAe = sum(JBTe.*mask,3);
FRIC = sum(JFT.*mask, 3);
WIND = sum(JFW.*mask, 3);

maskF = squeeze((RhoF > iso-delt) & (RhoF < iso+delt));
DIAs_F = sum(JBTs_FILT.*maskF, 3);
DIAe_F = sum(JBTe_FILT.*maskF, 3);
FRIC_F = sum(JFT_FILT.*maskF, 3);
WIND_F = sum(JFW_FILT.*maskF, 3);

T = squeeze(nanmean(Tf(:,:,end,:),4));

ml = double(T~=0);
ml(~ml) = NaN;
DIAs= DIAs.*ml;
DIA = DIA.*ml;
DIAe = DIAe.*ml;
FRIC = FRIC.*ml;
WIND = WIND.*ml;




%% LOAD SHORELINE
shorelines = gshhs('~/GSHH/gshhs_i.b');%_h for high, _c for crude, _i for intermediate, from https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
levels = [shorelines.Level];
land = (levels==1);
%%
xl = 600:1602;
yl = 1:550;
lats = lat(xl,yl); lons = lon(xl, yl);
lonf = [lons(1,:) lons(:,end).' fliplr(lons(end,:)) flipud(lons(:,1)).'];
latf = [lats(1,:) lats(:,end).' fliplr(lats(end,:)) flipud(lats(:,1)).'];

gshow = true;
cl = 5e-12;
cl = 1e-11; % For no hkpp bl mask
mc = .9; %coast color
gap = [.05 .01], margh = .1; margw=.1;
figure
subtightplot(2,2,1, gap, margh, margw);
axesm('mercator','MapLatLimit', [26.75 47], 'MapLonLimit', [-85 -56.5], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5,'MLabelParallel', 'South',...
    'PLineLocation', 5, 'MLineLocation', 5, 'Frame', 'on');
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(lat(:,:)),squeeze(lon(:,:)),(squeeze(T(:,:,end,1)))); shading interp
contourm(lat, lon, squeeze(T(:,:,end,1)), 0:1:30, 'k')
hold on
hf = fillm(latf, lonf, 'k');
set(hf, 'FaceColor', 'none');
set(hf, 'LineWidth', 1, 'LineStyle', '--');
hold off
if gshow; geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc); end
set(gca, 'clim', [10 28]);
t= textm(45, 360-82.5,0, '$T$', 'FontSize', 20, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');

subtightplot(2,2,2, gap, margh, margw);
axesm('mercator','MapLatLimit', [26.75 47], 'MapLonLimit', [-85 -56.5], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5,'MLabelParallel', 'South',...
    'PLineLocation', 5, 'MLineLocation', 5, 'Frame', 'on');
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(lat(:,:)),squeeze(lon(:,:)),(squeeze(-DIAs(:,:,end,1)))); shading interp
contourm(lat, lon, squeeze(T(:,:,end,1)), 0:1:30, 'k')
hold off
if gshow; geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc); end
% geoshow(map,refvec, 'DisplayType', 'texturemap')
set(gca, 'clim', [-1 1].*cl);
t= textm(45, 360-80,0, '$-J_D^{SURF}$', 'FontSize', 20, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');

subtightplot(2,2,3, gap, margh, margw);
axesm('mercator','MapLatLimit', [26.75 47], 'MapLonLimit', [-85 -56.5], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5,'MLabelParallel', 'South',...
    'PLineLocation', 5, 'MLineLocation', 5, 'Frame', 'on');
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(lat(:,:)),squeeze(lon(:,:)),squeeze(-FRIC-DIAe)); shading interp
contourm(lat, lon, squeeze(T(:,:,end,1)), 0:1:30, 'k')
hold off
if gshow; geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc); end
set(gca, 'clim', [-1 1].*cl);
t= textm(45, 360-80,0, '$-J^{SM}$', 'FontSize', 20, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');

subtightplot(2,2,4, gap, margh, margw);

axesm('mercator','MapLatLimit', [26.75 47], 'MapLonLimit', [-85 -56.5], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5,'MLabelParallel', 'South',...
    'PLineLocation', 5, 'MLineLocation', 5, 'Frame', 'on');
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(lat(:,:)),squeeze(lon(:,:)),squeeze(-WIND)); shading interp
contourm(lat, lon, squeeze(T(:,:,end,1)), 0:1:30, 'k')
hold off
if gshow; geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc); end
set(gca, 'clim', [-1 1].*cl);
t= textm(45, 360-80,0, '$-J_F^{WIND}$', 'FontSize', 20, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');

cb = colorbar;
ylabel(cb, 'PV Flux $[m s^{-4}]$', 'FontSize', 20, 'Interpreter', 'Latex');
set(gca, 'FontSize', 16)

colormap(cptcmap('BlueWhiteOrangeRed.cpt'));

set(gcf, 'Color','w', 'Position', [  307         159        1043         781]);

set(cb, 'Position', [0.9067    0.3093    0.0253    0.3752]);
pos = get(cb, 'Position');
haxes = axes('position', pos, 'color', 'none', 'ylim', [10 28], 'xtick', [], 'YTick', [10 19 28]);
ylabel('T $[^{\circ}C]$', 'FontSize', 16, 'Interpreter', 'Latex');

%% DECIMATED
reducF = squeeze(nansum(nansum((FRIC_F+DIAe_F).*1./pm.*1./pn)))./squeeze(nansum(nansum((FRIC+DIAe).*1./pm.*1./pn)));
reducD = squeeze(nansum(nansum((DIAs_F).*1./pm.*1./pn)))./squeeze(nansum(nansum((DIAs).*1./pm.*1./pn)));

xl = 600:1602;
yl = 1:550;
lats = lat(xl,yl); lons = lon(xl, yl);
lonf = [lons(1,:) lons(:,end).' fliplr(lons(end,:)) flipud(lons(:,1)).'];
latf = [lats(1,:) lats(:,end).' fliplr(lats(end,:)) flipud(lats(:,1)).'];

gshow = true;
cl = 5e-12;
cl = 1e-11; % For no hkpp bl mask
mc = .9; %coast color
gap = [.05 .01], margh = .1; margw=.1;
figure
subtightplot(2,2,1, gap, margh, margw);
axesm('mercator','MapLatLimit', [26.75 47], 'MapLonLimit', [-85 -56.5], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5,'MLabelParallel', 'South',...
    'PLineLocation', 5, 'MLineLocation', 5, 'Frame', 'on');
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(lat(:,:)),squeeze(lon(:,:)),(squeeze(-DIAs(:,:,end,1)))); shading interp
contourm(lat, lon, squeeze(T(:,:,end,1)), 0:1:30, 'k')
hold off
if gshow; geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc); end
% geoshow(map,refvec, 'DisplayType', 'texturemap')
set(gca, 'clim', [-1 1].*cl);
t= textm(44, 360-80,0, {'$-J_D^{SURF}$','$1.5$ km'}, 'FontSize', 20, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');


subtightplot(2,2,2, gap, margh, margw);
axesm('mercator','MapLatLimit', [26.75 47], 'MapLonLimit', [-85 -56.5], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5,'MLabelParallel', 'South',...
    'PLineLocation', 5, 'MLineLocation', 5, 'Frame', 'on');
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(lat(:,:)),squeeze(lon(:,:)),(squeeze(-(FRIC+DIAe)))); shading interp
contourm(lat, lon, squeeze(T(:,:,end,1)), 0:1:30, 'k')
hold off
if gshow; geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc); end
% geoshow(map,refvec, 'DisplayType', 'texturemap')
set(gca, 'clim', [-1 1].*cl);
t= textm(44, 360-80,0,{ '$-J^{SM}$';'$1.5$ km'}, 'FontSize', 20, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');

subtightplot(2,2,3, gap, margh, margw);
axesm('mercator','MapLatLimit', [26.75 47], 'MapLonLimit', [-85 -56.5], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5,'MLabelParallel', 'South',...
    'PLineLocation', 5, 'MLineLocation', 5, 'Frame', 'on');
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(lat(:,:)),squeeze(lon(:,:)),squeeze(-DIAs_F)); shading interp
% contourm(lat, lon, squeeze(nanmean(Tfilt(:,:,:),3)), 0:1:30, 'k')
contourm(lat, lon, squeeze(T(:,:,end,1)), 0:1:30, 'k')

hold off
if gshow; geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc); end
set(gca, 'clim', [-1 1].*cl);
t= textm(44, 360-80,0, {'$-J_D^{SURF}$'; '$10.5$ km'}, 'FontSize', 20, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');
t = textm(29, 360-62.5, 0, ['$\frac{LR}{HR} = ', num2str(reducD, 2), '$'], 'FontSize', 18, 'Color', 'k',...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');

subtightplot(2,2,4, gap, margh, margw);

axesm('mercator','MapLatLimit', [26.75 47], 'MapLonLimit', [-85 -56.5], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5,'MLabelParallel', 'South',...
    'PLineLocation', 5, 'MLineLocation', 5, 'Frame', 'on');
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(lat(:,:)),squeeze(lon(:,:)),squeeze(-(FRIC_F+DIAe_F))); shading interp
% contourm(lat, lon, squeeze(nanmean(Tfilt(:,:,:),3)), 0:1:30, 'k')
contourm(lat, lon, squeeze(T(:,:,end,1)), 0:1:30, 'k')

hold off
if gshow; geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc); end
set(gca, 'clim', [-1 1].*cl);
t= textm(44, 360-80,0, {'$-J^{SM}$';'$10.5$ km'}, 'FontSize', 20, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');
t = textm(29, 360-62.5, 0, ['$\frac{LR}{HR} = ', num2str(reducF, 2), '$'], 'FontSize', 18, 'Color', 'k',...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');

cb = colorbar;
ylabel(cb, 'PV Flux $[m s^{-4}]$', 'FontSize', 20, 'Interpreter', 'Latex');
set(gca, 'FontSize', 16)

colormap(cptcmap('BlueWhiteOrangeRed.cpt'));

set(gcf, 'Color','w', 'Position', [  307         159        1043         781]);

set(cb, 'Position', [0.905   0.3093    0.0253    0.3752]);
pos = get(cb, 'Position');
% haxes = axes('position', pos, 'color', 'none', 'ylim', [10 28], 'xtick', [], 'YTick', [10 19 28]);
% ylabel('T $[^{\circ}C]$', 'FontSize', 16, 'Interpreter', 'Latex')