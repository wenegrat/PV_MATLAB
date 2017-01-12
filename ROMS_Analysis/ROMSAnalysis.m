% NOTING all things to be looked at later:
% XXX - grid is spherical.
% XXX - Probably should directly calculate Pressure gradients (not press)
pardir = '/data/thomas/jacob13/GulfStream/NESEA/';
basepath = [pardir 'HIS/'];

files = dir([basepath,'*.nc']);
zl = 16:50;
xl = 1100:1400;
yl = 150:400;
slice =  {[xl(1) xl(end)], [yl(1) yl(end)], [zl(1) zl(end)], 0};

% Load constants/Initial Params
path = [basepath, files(1).name];
time = ncread(path, 'ocean_time');
nt = length(time);
ts = time(2)-time(1);


f = ncread('/data/thomas/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc', 'f');
f = f(xl, yl);
rho0 = 1027.4;
g = 9.81;

% FULL DOMAIN

T = GetVarROMS(path, 0, {'temp', '(1)'}, {0, 0,[48 50], 0});
 
pcolor(squeeze(T(:,:,end,1)).'); shading interp
hold on
rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
hold off

pm = ncread(path, 'pm');
pm = pm(xl, yl);
% pmmetric = repmat(pm, [1 1 nz nt]);

pn = ncread(path, 'pn');
pn = pn(xl,yl);
% pnmetric = repmat(pn, [1 1 nz nt]);

%%
% [nx, ny] = size(pm);
% PV = NaN(nx, ny, nz, 200);
% Tf = PV;
% 
% counter = 1;
% for i=1:length(files)-2
%     disp([num2str(i), '/', num2str(length(files)-1)]);
%     name = files(i).name;
%     path = [basepath, name];
%     
%     U = GetVarROMS(path, 0, {'u', '(1)'}, slice);
%     V = GetVarROMS(path, 0, {'v', '(1)'}, slice);
%     Vx = DrvROMS(pm, V, 'x');
%     Uy = DrvROMS(pn, U, 'y');
%     ZETA = Vx - Uy;
%     
%     T = GetVarROMS(path, 0, {'temp', '(1)'}, slice);
%     Tf(:,:,:,counter:counter+1) = T;
%     
%     S = GetVarROMS(path, 0, {'salt', '(1)'}, slice);
%     rho = roms_eos(T, S, 0); % XXX - referencing to 0 db.
%     b = -g*rho./rho0;
%     bz = DrvROMS(zmetric, b, 'z');
%     bx = DrvROMS(pm, b, 'x');
%     by = DrvROMS(pn, b, 'y');
%     
%     uz = DrvROMS(zmetric, U, 'z');
%     vz = DrvROMS(zmetric, V, 'z');
%     PV(:,:,:,counter:counter+1) = -vz.*bx + uz.*by+(repmat(f, [1 1 nz nt])+ZETA).*bz;
%     counter = counter + 2;
% end

%% Calculate Flux Terms
%need du/dt + ADV + CORI + PRESS
nt = 30;
[nx ny] = size(pm);
UPrev = NaN;
VPrev = NaN;
BPrev = NaN;
DuDt = NaN(nx, ny, nz, nt);
DvDt = DuDt;
UADV = DuDt;
VADV = DuDt;
UCOR = DuDt;
VCOR = DuDt;
DbDt = DuDt;
UPRE = DuDt;
VPRE = DuDt;
JFz = DuDt;
JFx = DuDt;
JFy = DuDt;
JBz = DuDt;
JBx = DuDt;
JBy = DuDt;
JAz = DuDt;
JAx = DuDt;
JAy = DuDt;
BADV = DuDt;
Q = DuDt;
Bx = DuDt;
By = DuDt;
Bz = DuDt;
zw = NaN(nx, ny, nz+1, nt);
OMEGAX = DuDt;
OMEGAY = DuDt;
OMEGAZ = DuDt;

    name = files(1).name;
    path = [basepath, name];
    
theta_s = ncreadatt(path, '/', 'theta_s');
theta_b = ncreadatt(path, '/', 'theta_b');
hc = ncreadatt(path, '/', 'hc');

h = ncread(path, 'h');
% x = ncread(path, 'lon_rho');
% y = ncread(path, 'lat_rho');
h = h(xl, yl);
% x = x(xl, yl); y = y(xl, yl);
for i=1:nt;((length(files)-2)*2);

    disp([num2str(i), '/', num2str((length(files)-2)*2)]);

    fileind = ceil(i/2);
    name = files(fileind).name;
    path = [basepath, name];
    
    Eta = GetVarROMS(path, 0, {'zeta', '(1)'}, sliceT);

    z = compZ(path, 0, Eta, theta_s, theta_b, hc, h);
    z = z(:,:,zl);
 
    zwt = compZ(path, 1, Eta,  theta_s, theta_b, hc, h);
    zw(:,:,:,i) = zwt(:,:,zl(1):zl(end)+1,:);
    
%         zw = compZ(path, 1);
%         zw = zw(xl, yl, zl(1):zl(end)+1);
        [~, ~, nz] = size(z);
        zm = squeeze(nanmean(nanmean(z)));
%         zmetric = repmat(z, [1 1 1 nt]);

    
 

    
    sliceind = 2 - mod(i,2);
    sliceT = {slice{1}, slice{2}, slice{3},[sliceind sliceind]};
    U = GetVarROMS(path, 0, {'u', '(1)'}, sliceT);
    %Time
    DuDt(:,:,:,i) = (U - UPrev)./ts;
%     if (i>1)
%       DuDtA = (DuDt(:,:,:,i) + DuDt(:,:,:,i-1))./2; 
%     else
%         DuDtA = DuDt(:,:,:,i);
%     end
%     
    UPrev = U;
    %Adv
    Ux = DrvS(pm,z, U, 'x');
    Uy = DrvS(pn,z, U, 'y');
    Uz = DrvROMS(z, U, 'z');
    
    % Meridional Velocity
    V = GetVarROMS(path, 0, {'v', '(1)'}, sliceT);
    DvDt(:,:,:,i) = (V - VPrev)./ts;
%         if (i>1)
%       DvDtA = (DvDt(:,:,:,i) + DvDt(:,:,:,i-1))./2; 
%     else
%         DvDtA = DvDt(:,:,:,i);
%     end
    VPrev = V;
    %Adv V
    Vx = DrvS(pm,z, V, 'x');
    Vy = DrvS(pn,z, V, 'y');
    Vz = DrvROMS(z, V, 'z');
    
    Wz = -Ux - Vy;
    
    %Buoyancy
    T = GetVarROMS(path, 0, {'temp', '(1)'}, sliceT);
    S = GetVarROMS(path, 0, {'salt', '(1)'}, sliceT);
%     rho = roms_eos(T, S, 0); % XXX - referencing to 0 db.
%     rho = rho + rho0;
    
    rho = rho_eos(T, S, z);
    B = -g*rho./rho0;
    
    ztemp = permute(z, [3 1 2]);
    zwtemp = permute(squeeze(zw(:,:,1:end-1,i)), [3 1 2]);
    Ttemp = permute(T, [3 1 2]);
    Stemp = permute(S, [3 1 2]);
    [rhoa, bvf] = rho_eos(Ttemp,Stemp,ztemp, zwtemp, g, rho0);
    bvf = permute(bvf, [2 3 1]);
    DbDt(:,:,:,i) = (B - BPrev)./ts;
%     if (i>1)
%       DbDtA = (DbDt(:,:,:,i) + DbDt(:,:,:,i-1))./2; 
%     else
%         DbDtA = DbDt(:,:,:,i);
%     end
%     
    BPrev = B;
    Bx(:,:,:,i) = DrvS(pm, z, B, 'x');
    By(:,:,:,i) = DrvS(pn, z, B, 'y');
    % - XXX Should be adiabatically leveled first....
    Bz(:,:,:,i) = DrvROMS(z, B, 'z');
    Bz(:,:,:,i) = bvf;
    % Iterate for vertical integration terms
    
    omega = GetVarROMS(path, 0, {'omega', '(1)'}, sliceT);
    Zx = DrvROMS(pm, z, 'x'); % XXX-Is this the right type of derivative to take?
    Zy = DrvROMS(pn, z, 'y');
    W = U.*Zx + V.*Zy + omega;% see: https://www.myroms.org/forum/viewtopic.php?f=19&t=2139
    
%     W = NaN(nx, ny, nz);
    P = NaN(nx, ny, nz);
    for x =1:nx
        for y =1:ny
            mask = isfinite(Wz(x,y,:));
%             W(x, y, mask) = cumtrapz(squeeze(z(x,y,mask)), squeeze(Wz(x,y,mask)));
            zt = squeeze(z(x,y,:));
            bt = squeeze(B(x,y,:));
            mask = isfinite(bt);
            P(x, y, mask) = -rho0.*B(x, y, end).*Eta(x,y) + rho0.*flipud(cumtrapz(flipud(zt(mask)), flipud(bt(mask))));
        end
    end
    
    Px = DrvS(pm, z, P, 'x');
    Py = DrvS(pn, z, P, 'y');
    
    UADV(:,:,:,i) = U.*Ux + V.*Uy + W.*Uz;
    VADV(:,:,:,i) = U.*Vx + V.*Vy + W.*Vz;
    UCOR(:,:,:,i) = -repmat(f, [1 1 nz]).*V;
    VCOR(:,:,:,i) =  repmat(f, [1 1 nz]).*U;
    UPRE(:,:,:,i) = 1/rho0.*Px;
    VPRE(:,:,:,i) = 1/rho0.*Py;
    
    BADV(:,:,:,i) = U.*Bx(:,:,:,i) + V.*By(:,:,:,i) + W.*Bz(:,:,:,i);
    OMEGAX(:,:,:,i) = - Vz;
    OMEGAY(:,:,:,i) = Uz;
    OMEGAZ(:,:,:,i) = repmat(f, [1 1 nz]) + Vx - Uy;
    
    Q(:,:,:,i) = OMEGAX(:,:,:,i).*Bx(:,:,:,i) + OMEGAY(:,:,:,i).*By(:,:,:,i)+OMEGAZ(:,:,:,i).*Bz(:,:,:,i);
    
    JAz(:,:,:,i) = W.*Q(:,:,:,i); JAx(:,:,:,i) = U.*Q(:,:,:,i); JAy(:,:,:,i) = V.*Q(:,:,:,i);

end

%Shift time derivatives 
% XXX-Technically should also shift verticaly to account for z shifting...
DuDt(:,:,:,1:end-1) = (DuDt(:,:,:,1:end-1) + DuDt(:,:,:,2:end))./2;
DvDt(:,:,:,1:end-1) = (DvDt(:,:,:,1:end-1) + DvDt(:,:,:,2:end))./2;
DbDt(:,:,:,1:end-1) = (DbDt(:,:,:,1:end-1) + DbDt(:,:,:,2:end))./2;

for i=1:nt;((length(files)-2)*2);
    FTERMX = DuDt(:,:,:,i) + UADV(:,:,:,i) + UPRE(:,:,:,i)+UCOR(:,:,:,i);
    FTERMY = DvDt(:,:,:,i) + VADV(:,:,:,i) + VPRE(:,:,:,i)+VCOR(:,:,:,i);
    JFz(:,:,:,i) = Bx(:,:,:,i).*FTERMY - By(:,:,:,i).*FTERMX;
    JFx(:,:,:,i) = -Bz(:,:,:,i).*FTERMY; 
    JFy(:,:,:,i) = Bz(:,:,:,i).*FTERMX;
    
    D = DbDt(:,:,:,i) + BADV(:,:,:,i);

    JBx(:,:,:,i) = -OMEGAX(:,:,:,i).*D; JBy(:,:,:,i) = -OMEGAY(:,:,:,i).*D; JBz(:,:,:,i) = -OMEGAZ(:,:,:,i).*D;
    
end


%%
mask = ones(size(Q));
dz = diff(zw, 1, 3);
gridvol = abs(repmat(repmat(1./pm, [1 1 nz]).*repmat(1./pn, [1 1 nz]), [1 1 1 nt])).*abs(dz);
vol = squeeze(sum(sum(sum(mask.*gridvol))));

IntegrateQTerms;

%Vertical Terms
dx = repmat(1./pm, [1 1 nt]);
dy = repmat(1./pn, [1 1 nt]);
[JFz_t, ~    ] = areaIntegrateJVecsROMS(squeeze(JFz(:,:,end-1,:)), squeeze(mask(:,:,end-1,:)), dx, dy, ts, vol);
[JFz_b, ~    ] = areaIntegrateJVecsROMS(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx, dy, ts, vol);
JFzA = JFz_t - JFz_b;
[JBz_t, ~    ] = areaIntegrateJVecsROMS(squeeze(JBz(:,:,end-1,:)), squeeze(mask(:,:,end-1,:)), dx, dy, ts, vol);
[JBz_b, ~    ] = areaIntegrateJVecsROMS(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx, dy, ts, vol);
JBzA = JBz_t - JBz_b;
[JAzA, ~] = areaIntegrateJVecsROMS(squeeze(JAz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx, dy, ts, vol);
JAzA = -JAzA;

%Zonal Terms
dyl = repmat(squeeze(1./pn(3,:).'), [1 nz nt]);
dyr = repmat(squeeze(1./pn(end-2,:).'), [1 nz nt]);
% dz = diff(zw,1,3);
% dzl = repmat(squeeze(dz(3,:,:)), [1 1 nt]);
% dzr = repmat(squeeze(dz(end-2,:,:)), [1 1 nt]);
dzl = squeeze(dz(3,:,:,:));
dzr = squeeze(dz(end-2,:,:,:));
[JFx_l, ~    ] = areaIntegrateJVecsROMS(squeeze(JFx(3,:,:,:)), squeeze(mask(3,:,:,:)),  dyl, dzl, ts, vol);
[JFx_r, ~    ] = areaIntegrateJVecsROMS(squeeze(JFx(end-2,:,:,:)), squeeze(mask(end-2,:,:,:)), dyr, dzr, ts, vol);
JFxA = JFx_r - JFx_l;
[JBx_l, ~    ] = areaIntegrateJVecsROMS(squeeze(JBx(3,:,:,:)), squeeze(mask(3,:,:,:)),  dyl, dzl, ts, vol);
[JBx_r, ~    ] = areaIntegrateJVecsROMS(squeeze(JBx(end-2,:,:,:)), squeeze(mask(end-2,:,:,:)), dyr, dzr, ts, vol);
JBxA = JBx_r - JBx_l;
[JAx_l, ~    ] = areaIntegrateJVecsROMS(squeeze(JAx(3,:,:,:)), squeeze(mask(3,:,:,:)),  dyl, dzl, ts, vol);
[JAx_r, ~    ] = areaIntegrateJVecsROMS(squeeze(JAx(end-2,:,:,:)), squeeze(mask(end-2,:,:,:)), dyr, dzr, ts, vol);
JAxA = JAx_r - JAx_l;

%Meridional Terms
dxf = repmat(squeeze(1./pm(:,3)), [1 nz nt]);
dxb = repmat(squeeze(1./pm(:,end-2)), [1 nz nt]);
% dzf = repmat(squeeze(dz(:,3,:)), [1 1 nt]);
% dzb = repmat(squeeze(dz(:,end-2,:)), [1 1 nt]);
dzf = squeeze(dz(:,3,:,:));
dzb = squeeze(dz(:,end-2,:,:));
[JFy_f, ~    ] = areaIntegrateJVecsROMS(squeeze(JFy(:,3,:,:)), squeeze(mask(:,3,:,:)),  dxf, dzf, ts, vol);
[JFy_b, ~    ] = areaIntegrateJVecsROMS(squeeze(JFy(:,end-2,:,:)), squeeze(mask(:,end-2,:,:)), dxb, dzb, ts, vol);
JFyA = JFy_b - JFy_f;
[JBy_f, ~    ] = areaIntegrateJVecsROMS(squeeze(JBy(:,3,:,:)), squeeze(mask(:,3,:,:)),  dxf, dzf, ts, vol);
[JBy_b, ~    ] = areaIntegrateJVecsROMS(squeeze(JBy(:,end-2,:,:)), squeeze(mask(:,end-2,:,:)), dxb, dzb, ts, vol);
JByA = JBy_b - JBy_f;
[JAy_f, ~    ] = areaIntegrateJVecsROMS(squeeze(JAy(:,3,:,:)), squeeze(mask(:,3,:,:)),  dxf, dzf, ts, vol);
[JAy_b, ~    ] = areaIntegrateJVecsROMS(squeeze(JAy(:,end-2,:,:)), squeeze(mask(:,end-2,:,:)), dxb, dzb, ts, vol);
JAyA = JAy_b - JAy_f;
%%
figure
plot(Qa, 'LineWidth', 2);
hold on
plot(-JFzA);
plot(-JBzA);
plot(-JFxA - JFyA);
plot(-JBxA - JByA);
plot(-JAxA - JAyA-JAzA);
plot(-(JFzA + JBzA+JFxA+JFyA+JBxA+JByA+JAxA+JAyA+JAzA), '--', 'LineWidth', 2);
plot(-JFzA-JBzA);
hold off
legend('\Delta Q','-JFz', '-JBz', '-JFh', '-JBh', '-JA', 'Sum', 'Location', 'SouthWest')
grid on
xlabel('Timestep (12 hours)');
set(gca, 'FontSize', 18);

%%
forcepath = [pardir 'nesea_frc.nc'];
Qor = ncread(forcepath,'shflux'); % Monthly Values W/m^2
EPr = ncread(forcepath,'swflux').*100./86400; %Freshwater flux m/s
tau_u = ncread(forcepath,'sustr'); %N/m^2;
tau_v = ncread(forcepath, 'svstr'); %N/m^2

climtime = datenum(2012,1:12, 1);
modeltime = datenum(2012,1,1, 1:12:(nt*12),0,0);

Qo = NaN(nx, ny, nt);
EP = Qo;

for x = 1:nx;
    for y = 1:ny
Qo(x, y,:) = interp1(climtime, squeeze(Qor(x,y,:)), modeltime, 'pchip');
EP(x,y,:) = interp1(climtime, squeeze(EPr(x,y,:)), modeltime, 'pchip');
    end
end

% alpha = 

%%
plot(-JAzA);
hold on
plot(-JAxA);
plot(-JAyA);
hold off