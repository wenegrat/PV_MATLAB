% NOTING all things to be looked at later:
% XXX - grid is spherical.
% XXX - Probably should directly calculate Pressure gradients (not press)

%PARAMETERS TO CHANGE'
% pardir = '/groups/thomas1/jacob13/GulfStream/NESEA/';
% basepath = [pardir 'HIS/'];
% path1 = '/groups/thomas1/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc';
% ntp = 2;
% nt = 50;
% forcepath = [pardir 'nesea_frc.nc'];
% modeltime = datenum(2012,1,1, 1:12:(nt*12),0,0);

pardir = '/data/thomas/jacob13/GULFZ/';
basepath = [pardir 'HIS/'];
path1 = [pardir 'gulfz_grd.nc'];
forcepath = [pardir 'gulfz_frc.nc'];
ntp = 5;
nt = 345;
modeltime = datenum(2012,8,(26++offset):1:(26+offset+(nt)-1), 0,0,0);

%
files = dir([basepath,'*.nc']);
zl = 45:50;
nz = length(zl);
xl = 1:1602;
yl = 1:922;
slice =  {[xl(1) xl(end)], [yl(1) yl(end)], [zl(1) zl(end)], 0};

% Load constants/Initial Params
path = [basepath, files(1).name];
time = ncread(path, 'ocean_time');
% nt = length(time);
ts = time(2)-time(1);

path = path1;
f = ncread(path, 'f');
f = f(xl, yl);
rho0 = 1027.4;
g = 9.81;
pm = ncread(path, 'pm');
pm = pm(xl, yl);
% pmmetric = repmat(pm, [1 1 nz nt]);

pn = ncread(path, 'pn');
pn = pn(xl,yl);

h = ncread(path, 'h');
subplot(2,1,2)
pcolor(h.'); shading interp
hold on
rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
hold off
colorbar

% FULL DOMAIN
path = [basepath, files(50).name];

T = GetVarROMS(path, 0, {'temp', '(1)'}, {0, 0,0, 0});

subplot(2,1,1)
pcolor(squeeze(T(:,:,24,1)).'); shading interp
hold on
rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
hold off
colorbar;




%% Calculate Flux Terms
%need du/dt + ADV + CORI + PRESS
offset = 0; % Sets the 'start date'
[nx ny] = size(pm);

DuDt = NaN(nx, ny, nz, nt);

Bx = DuDt;
By = DuDt;

Tf = DuDt;
Sf = Tf;
Rho = DuDt;
hkpp = NaN(nx, ny, nt);


    name = files(1).name;
    path = [basepath, name];
    
theta_s = ncreadatt(path, '/', 'theta_s');
theta_b = ncreadatt(path, '/', 'theta_b');
hc = ncreadatt(path, '/', 'hc');

% h = ncread(path, 'h');
% x = ncread(path, 'lon_rho');
% y = ncread(path, 'lat_rho');
h = h(xl, yl);
% x = x(xl, yl); y = y(xl, yl);
for i=1:nt;
    % 1) Housekeeping
    disp([num2str(i), '/', num2str((length(files)-ntp)*ntp)]);
    fileind = ceil((i+offset)/ntp);
    name = files(fileind).name;
    path = [basepath, name];
    sliceind = ntp - mod(i,ntp);
    sliceT = {slice{1}, slice{2}, slice{3},[sliceind sliceind]};
    
    % 2) Calculate z coordinates at this timestep.
    Eta = GetVarROMS(path, 0, {'zeta', '(1)'}, sliceT);

    hkpp(:,:,i) = GetVarROMS(path, 0, {'hbls', '(1)'}, sliceT);
    
    z = compZ(path, 0, Eta, theta_s, theta_b, hc, h);
    z = z(:,:,zl);
    zwt = compZ(path, 1, Eta,  theta_s, theta_b, hc, h);
%     zw(:,:,:,i) = zwt(:,:,zl(1):zl(end)+1,:);  % Save zw for later dz calc
%         zw = compZ(path, 1);
%         zw = zw(xl, yl, zl(1):zl(end)+1);
    [~, ~, nz] = size(z);
    zm = squeeze(nanmean(nanmean(z)));
    % Omega is the vertical velocity in s-coordinates


    
    % 5) Load Buoyancy Terms
    T = GetVarROMS(path, 0, {'temp', '(1)'}, sliceT);
     Tf(:,:,:,i) = T;
    S = GetVarROMS(path, 0, {'salt', '(1)'}, sliceT);
     Sf(:,:,:,i) = S;
%     rho = roms_eos(T, S, 0); % XXX - referencing to 0 db.
%     rho = rho + rho0;
    rho = rho_eos(T, S, 0); % CROCO function (checked for consistency 1/12/17)
    Rho(:,:,:,i) = rho;
    B = -g*rho./rho0; % In-Situ B

   

    Bx(:,:,:,i) = DrvS(pm, z, B, 'x');
    By(:,:,:,i) = DrvS(pn, z, B, 'y');
    % - XXX Should be adiabatically leveled first....
%     Bz(:,:,:,i) = bvf;
    % Iterate for vertical integration terms
end


%%

Qor = ncread(forcepath,'shflux'); % Monthly Values W/m^2
EPr = ncread(forcepath,'swflux')./(100.*86400); %Freshwater flux m/s
tau_um = ncread(forcepath,'sustr'); %N/m^2;
tau_vm = ncread(forcepath, 'svstr'); %N/m^2
Qor = Qor(xl, yl,:);
EPr = EPr(xl, yl,:);
tau_u(xl(2:end),yl, :) = tau_um;%(xl, yl,:);
tau_v(xl, yl(2:end),:) = tau_vm;%(xl, yl(2:end),:);
SST = ncread(forcepath, 'SST');
SST = SST(xl, yl,:);
Qor = repmat(Qor, [1 1 2]);
EPr = repmat(EPr, [1 1 2]);
tau_u = repmat(tau_u, [1 1 2]);
tau_v = repmat(tau_v, [1 1 2]);
SST = repmat(SST, [1 1 2]);
climtime = datenum(2012,1, 15:30:(360*2));

% disp('here')
Qo = NaN(nx, ny, nt);
EP = Qo;
tx = Qo;
sst = Qo;
ty = Qo;
for x = 1:nx;
    for y = 1:ny
Qo(x, y,:) = interp1(climtime, squeeze(Qor(x,y,:)), modeltime, 'pchip');
EP(x,y,:) = interp1(climtime, squeeze(EPr(x,y,:)), modeltime, 'pchip');
tx(x,y,:) = interp1(climtime, squeeze(tau_u(x,y,:)), modeltime, 'pchip');
ty(x,y,:) = interp1(climtime, squeeze(tau_v(x,y,:)), modeltime, 'pchip');
sst(x,y,:) = interp1(climtime, squeeze(SST(x,y,:)), modeltime, 'pchip');
    end
end

Tmean = squeeze(nanmean(nanmean(nanmean(Tf(:,:,end,:), 4))));
Smean = squeeze(nanmean(nanmean(nanmean(Sf(:,:,end,:), 4))));

%% GEN THEORY SCALINGS
del = 1;

rhot = rho_eos(squeeze(Tf(:,:,end,:)), squeeze(Sf(:,:,end,:)), 0);
alpha = 1./rhot.*(rho_eos(squeeze(Tf(:,:,end,:))+del, squeeze(Sf(:,:,end,:)), 0)-rho_eos(squeeze(Tf(:,:,end,:))-del, squeeze(Sf(:,:,end,:)), 0))./(2*del);
beta = 1./rhot.*(rho_eos(squeeze(Tf(:,:,end,:)), squeeze(Sf(:,:,end,:))+del, 0)-rho_eos(squeeze(Tf(:,:,end,:)), squeeze(Sf(:,:,end,:))-del, 0))./(2*del);

Cp = 3994; % XXX - Get ROMS Value
ssta =  squeeze(Tf(:,:,end,:))-sst;
Bo = g*alpha.*(Qo-30*ssta)./(rho0*Cp) +  g*beta.*EP.*squeeze(Sf(:,:,end,:));
% Bo = g*alpha.*(Qo)./(rho0*Cp) +  g*beta.*EP.*squeeze(Sf(:,:,end,:));
% zf = repmat(zw, [1 1 1 nt]);
zf = zw(:,:,1:end-1,:);
H = hkpp; H(:,:,1) = H(:,:,2);
% H(mldepth>H) = mldepth(mldepth>H);
% H = squeeze(mldepth);
% H(H<12) = 12;
% ms = double(zf > -permute(repmat(H, [1 1 1 nz]), [1 2 4 3]));
% ms(ms ==0) = NaN;
% M2 = squeeze(squeeze(nanmean(ms.*Bx(:,:,:,:), 3)).^2 + squeeze(nanmean(ms.*By(:,:,:,:), 3)).^2); %  XX- not averaging over BL...
% m = ~isfinite(M2);  
% Ma = squeeze(Bx(:,:,end-1,:)).^2 + squeeze(By(:,:,end-1,:)).^2;
% M2(m) = Ma(m); % Use uppermost gradient value if H is too small...

M2 = squeeze( nanmean(Bx(:,:,:,:),3).^2 + nanmean(By(:,:,:,:),3).^2);

tmag = abs(tx+1i*ty);
% de = sqrt(2*.1*H.*sqrt(tmag./rho0)./repmat(f,[1 1 nt]));

% H(H<de) = de(H<de);
H(H<6) = 6;
% H = repmat(nanmean(nanmean(H)), [nx ny 1]);
% JBT = repmat(f, [1 1 nt]).*Bo./H + 0.15*M2.*H;
% JBT(hkpp < 20) = NaN;
% JBTs = JBT - 0.15.*M2.*H;
JFT = -0.2*M2.*H;
% H(H<24) = 24;
Hm =H;
% Hm(Hm<14 & Bo<0) = 14; % Logic here is if the boundary layer is shallow and flux is into ML, correct scale is absorption scale.
JBTs = repmat(f, [1 1 nt]).*Bo./Hm;
% JBTs(JBTs<0) = 0;
JBTe = 0.15.*M2.*H;
JBT = JBTs + JBTe;

% disp('here')

R = 0.05.*(H.^2.*M2)./(Bo.*repmat(f, [1 1 nt]));

JFW = -tx./(rho0*H).*squeeze(By(:,:,end-2,:)) + ty./(rho0*H).*squeeze(Bx(:,:,end-2,:));



%%
iso = 1026;
delt = 0.2; 
mask = squeeze((Rho(:,:,end-1,:) > iso-delt) & (Rho(:,:,end-1,:) < iso+delt));
mask = mask &hkpp>100;
DIA = sum(JBT.*mask, 3);
DIAs = sum(JBTs.*mask,3);
DIAe = sum(JBTe.*mask,3);
FRIC = sum(JFT.*mask, 3);
WIND = sum(JFW.*mask, 3);
T = squeeze(nanmean(Tf(:,:,end,:),4));

ml = double(T~=0);
ml(~ml) = NaN;
DIAs= DIAs.*ml;
DIA = DIA.*ml;
DIAe = DIAe.*ml;
FRIC = FRIC.*ml;
WIND = WIND.*ml;
%%
cl = 5e-12;
figure
subplot(4,1,1);
pcolor(T.'); shading interp
colorbar;
hold on
contour(T.', [17 19], 'k');
hold off

subplot(4,1,2)
pcolor(-DIAs.'); shading interp
hold on
contour(T.', [17 19], 'k');
hold off
colorbar
set(gca, 'clim', [-1 1].*cl);

subplot(4,1,3)
pcolor(-(FRIC+DIAe).') ; shading interp
colorbar
hold on
contour(T.', [17 19], 'k');
hold off
set(gca, 'clim', [-1 1].*cl);

subplot(4,1,4)
pcolor(-WIND.'); shading interp
colorbar
hold on
contour(T.', [17 19], 'k');
hold off
set(gca, 'clim', [-1 1].*cl);

%%
lat = ncread(path1, 'lat_rho');
lon = ncread(path1, 'lon_rho');

m_proj('mercator', 'longitudes', [-90 0], 'latitudes', [15 45]);
m_coast('patch',[.7 .7 .7]);
m_grid('box', 'fancy', 'fontsize', 26);
% m_contourf([LON; LON-360], LAT, [squeeze(MT(:,:,1)).' squeeze(MT(:,:,1)).'], 10);

m_pcolor(lon(:,1), lat(1,:), -DIA.');
