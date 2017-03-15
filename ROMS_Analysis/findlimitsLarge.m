function [xl, yl, offset, duration] = findlimitsLarge(xls, yls, nt)

pathsm = '/groups/thomas1/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc';
latst = ncread(pathsm, 'lat_rho');
latSm = latst(xls, yls);
lonst = ncread(pathsm, 'lon_rho');
lonSm = lonst(xls, yls);

pardir = '/data/thomas/jacob13/GULFZ/';
basepath = [pardir 'HIS/'];
pathL = [pardir 'gulfz_grd.nc'];

latsL = ncread(pathL, 'lat_rho');
lonsL = ncread(pathL, 'lon_rho');

disp(['Lat Limit Sm: ', num2str(latSm(1,1)),'-', num2str(latSm(end,1))]);

dist = sqrt((latsL - latSm(1,1)).^2 + (lonsL - lonSm(1,1)).^2 );
[minVal, minIndex] = min(dist(:));

cbl = cell(1, ndims(dist));
[cbl{:}] = ind2sub(size(dist), minIndex);

dist = sqrt((latsL - latSm(end,end)).^2 + (lonsL - lonSm(end,end)).^2 );
[minVal, minIndex] = min(dist(:));

ctr = cell(1, ndims(dist));
[ctr{:}] = ind2sub(size(dist), minIndex);

xl = cbl{1}:ctr{1};
yl = cbl{2}:ctr{2};


%% TIME

files = dir(['/groups/thomas1/jacob13/GulfStream/NESEA/HIS/','*.nc']);
path = ['/groups/thomas1/jacob13/GulfStream/NESEA/HIS/', files(1).name];
tf = ncread(path, 'ocean_time');
time = ncread(path, 'ocean_time');
ts = time(2)-time(1);
modeltimeSm = tf(1):ts:((nt-1)*ts+tf(1));

files = dir(['/data/thomas/jacob13/GULFZ/HIS/','*.nc']);

tf = ncread([basepath, files(1).name], 'ocean_time');
path = [basepath, files(1).name];
time = ncread(path, 'ocean_time');
ts = time(2)-time(1);


offsets = 0:10000;

ind = find(tf(1) + offsets.*ts >= modeltimeSm(1), 1, 'first');
offset = offsets(ind);

ind = find(tf(1) + offsets.*ts <= modeltimeSm(end), 1, 'last');
duration = offsets(ind)-offset;
end