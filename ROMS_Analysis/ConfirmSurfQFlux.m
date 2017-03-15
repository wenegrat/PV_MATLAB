%%
% CAN REPRODUCE THE ISSUE IN NESEA using nt = 43; offset = 150; xl =
% 800:1050; yl = 300:500. Working theory is advection of GS meander into
% domain is not captured by fluxes at edge of isopycnal.

%PARAMETERS TO CHANGE'
pardir = '/groups/thomas1/jacob13/GulfStream/NESEA/';
basepath = [pardir 'HIS/'];
path1 = '/groups/thomas1/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc';
ntp = 2;
nt = 190;
nt = 400;
% nt = 210;
% nt = 43;
offset = 0;
% offset = 81; % use for NESEB comparison.
% nt = 35; % use for NESB comparison.

% offset = 0; % Use for NESEB
% nt = 418; % Use for NESEB

forcepath = [pardir 'nesea_frc.nc'];
bdom = false;

% % 
pardir = '/data/thomas/jacob13/GULFZ/';
basepath = [pardir 'HIS/'];
path1 = [pardir 'gulfz_grd.nc'];
forcepath = [pardir 'gulfz_frc.nc'];
ntp = 5;
bdom = true;
% % nt =  93;345;
% % offset =  126;



files = dir([basepath,'*.nc']);


% Load constants/Initial Params
path = [basepath, files(1).name];
time = ncread(path, 'ocean_time');
% nt = length(time);
ts = time(2)-time(1);


% startmodel = startdate + datenum(0, 0, 0, 0, 0, tf(1)-oti);
% dateoff = datenum(0, 0, offset);
% modeltime = inidate + dateoff + datenum(0, 0, 0:1:nt-1);


% modeltime = datenum(2012,8,(26+offset):1:(26+offset+(nt)-1), 0,0,0);


%

zl = 49:50;
nz = length(zl);
xl =  1:2002;100:600;1000:1150;1200:1300;
yl =  1:1602;100:500; 1:150;250:350;

% USE THIS FOR LARGE DOMAIN TO OVERLAP.
if bdom
[xl, yl, offset, nt] = findlimitsLarge(xl, yl, nt);
end

if offset>0
tf = ncread([basepath, files(ceil(offset./ntp)).name], 'ocean_time');
else
    tf = ncread([basepath, files(1).name], 'ocean_time');
end
modeltime = tf(1):ts:((nt-1)*ts+tf(1));

% xl = 1237:1336;
% yl = 295:395;
slice =  {[xl(1) xl(end)], [yl(1) yl(end)], [zl(1) zl(end)], 0};



path = path1;
f = ncread(path, 'f');
f = f(xl, yl);
rho0 = 1027.4;
g = 9.81;
Cp = 3994; % XXX - Get ROMS Value
zmin = 5.86;

% pmmetric = repmat(pm, [1 1 nz nt]);




latst = ncread(path1, 'lat_rho');
lats = latst(xl, yl);
lonst = ncread(path1, 'lon_rho');
lons = lonst(xl, yl);
disp(['Sub-domain borders (clockwise from bottom left): (', num2str(lats(1,1)), ',',num2str(lons(1,1)),')  (',...
    num2str(lats(1,end)),',',num2str(lons(1,end)),')  (', num2str(lats(end,end)),',',num2str(lons(end,end)),')  (',...
    num2str(lats(end,1)),',',num2str(lons(end,1)),')']);
% FULL DOMAIN
 
disp(['Starting DOY: ', num2str(tf(1)./86400)]);
disp(['Ending   DOY: ', num2str(modeltime(195)./86400)]);
disp(['TOTAL LENGTH: ', num2str((modeltime(195)-modeltime(1))./86400), ' days']);
%%
rho = NaN(length(xl), length(yl), nz, nt);
Tf = NaN(length(xl), length(yl), nt);
for i=1:nt
    disp([num2str(i), '/', num2str((length(files)-ntp)*ntp)]);
    fileind = ceil((i+offset)/ntp);
    name = files(fileind).name;
    path = [basepath, name];
    sliceind = mod(i,ntp);
    if sliceind==0; sliceind=ntp;end
    sliceT = {slice{1}, slice{2}, slice{3},[sliceind sliceind]};
    
    T = GetVarROMS(path, 0, {'temp', '(1)'}, sliceT);
%     Tf(:,:,:,i)=T; %% Might be worth keeping this later for plotting purposes...
    S = GetVarROMS(path, 0, {'salt', '(1)'}, sliceT);

    
    rho(:,:,:,i) = rho_eos(T, S, 0); % CROCO function (checked for consistency 1/12/17)
    Tf(:,:,i) = squeeze(T(:,:,end));

end
    mask = (1025.9<rho) & (1026.3>rho);
    mask = squeeze(mask(:,:,end,:));
    
%% LOAD FLUXES
Qor = ncread(forcepath,'shflux'); % Monthly Values W/m^2
Qor = Qor(xl, yl,:);
Qor = repmat(Qor, [1 1 2]);
SST = ncread(forcepath, 'SST');
SST = SST(xl, yl,:);
SST = repmat(SST, [1 1 2]);

% climtime = datenum(2012,1, 15:30:(360*2));
ds = 86400;
climtime = ((360+15)*ds):30*ds:(360*2*ds+360*ds);
% disp('here')
[nx, ny, ~] = size(Qor);

Qo = NaN(nx, ny, nt);


[X, Y, T] = ndgrid(xl, yl, climtime);
[Xm, Ym, mt] = ndgrid(xl, yl, modeltime);
QG = griddedInterpolant(X, Y, T, Qor, 'cubic', 'none');
SSTG = griddedInterpolant(X, Y, T, SST,'cubic', 'none');
sst = SSTG(Xm, Ym, mt);


Qo = QG(Xm, Ym, mt);


ssta = Tf - sst;
Qeff = Qo - 30*ssta;

% clear Qor EPr tau_u tau_v SST SWR


    
    %%
    
    figure
    plot(squeeze(nanmean(nanmean(Qo))));
    hold on
        plot(squeeze(min(min(Qo))), '--');
            plot(squeeze(nanmean(nanmean(Qeff))));

            plot(squeeze(nanmean(nanmean(Qeff(:,:,:).*mask))), ':');

    yt = get(gca, 'YTick');
    hold on
    plot(195./2.*ones(size(yt)), yt)
    hold off
    grid on
    legend('<Q_o>', 'min(Q_o)', '<Q_{Eff}>', '<Q_{Eff}>|_{MW}');