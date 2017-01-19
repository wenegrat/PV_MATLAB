% NOTING all things to be looked at later:
% XXX - grid is spherical.
% XXX - Probably should directly calculate Pressure gradients (not press)
pardir = '/data/thomas/jacob13/GulfStream/NESEA/';
basepath = [pardir 'HIS/'];

files = dir([basepath,'*.nc']);
zl = 15:50;
nz = length(zl);
xl = 950:1150;
yl = 150:400;
slice =  {[xl(1) xl(end)], [yl(1) yl(end)], [zl(1) zl(end)], 0};

% Load constants/Initial Params
path = [basepath, files(1).name];
time = ncread(path, 'ocean_time');
% nt = length(time);
ts = time(2)-time(1);


f = ncread('/data/thomas/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc', 'f');
f = f(xl, yl);
rho0 = 1027.4;
g = 9.81;

% FULL DOMAIN
path = [basepath, files(1).name];

T = GetVarROMS(path, 0, {'temp', '(1)'}, {0, 0,[48 50], 0});
figure
subplot(2,1,1)
pcolor(squeeze(T(:,:,end,1)).'); shading interp
hold on
rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
hold off
colorbar;

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

%%
nt= 190;
% pnmetric = repmat(pn, [1 1 nz nt]);
Ts = NaN(2002, 1602, nt); 
Hs = Ts;
for i=1:nt;length(files);
    
    disp([num2str(i), '/', num2str((length(files)-2)*2)]);
    fileind = ceil(i/2);
    name = files(fileind).name;
    path = [basepath, name];
    sliceind = 2 - mod(i,2);
    sliceT = {0, 0, [50 50],[sliceind sliceind]};
    Ts(:,:,i) = squeeze(GetVarROMS(path, 0, {'temp', '(1)'}, sliceT));
    Hs(:,:,i) = squeeze(GetVarROMS(path, 0, {'hbls', '(1)'}, sliceT));
end

%%
for i=1:nt
   pcolor(squeeze(Hs(:,:,i)).'); shading interp
   colorbar
   hold on
rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
hold off
title(num2str(i))
drawnow
 pause(0.001);
end
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
nt = 194;
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
Tf = DuDt;
Sf = DuDt;
Rho = DuDt;
zw = NaN(nx, ny, nz+1, nt);
hkpp = NaN(nx, ny, nt);
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
    % 1) Housekeeping
    disp([num2str(i), '/', num2str((length(files)-2)*2)]);
    fileind = ceil(i/2);
    name = files(fileind).name;
    path = [basepath, name];
    sliceind = 2 - mod(i,2);
    sliceT = {slice{1}, slice{2}, slice{3},[sliceind sliceind]};
    
    % 2) Calculate z coordinates at this timestep.
    Eta = GetVarROMS(path, 0, {'zeta', '(1)'}, sliceT);
    hkpp(:,:,i) = GetVarROMS(path, 0, {'hbls', '(1)'}, sliceT);
    
    z = compZ(path, 0, Eta, theta_s, theta_b, hc, h);
    z = z(:,:,zl);
    zwt = compZ(path, 1, Eta,  theta_s, theta_b, hc, h);
    zw(:,:,:,i) = zwt(:,:,zl(1):zl(end)+1,:);  % Save zw for later dz calc
%         zw = compZ(path, 1);
%         zw = zw(xl, yl, zl(1):zl(end)+1);
    [~, ~, nz] = size(z);
    zm = squeeze(nanmean(nanmean(z)));
    % Omega is the vertical velocity in s-coordinates
    omega = GetVarROMS(path, 0, {'omega', '(1)'}, sliceT);
    Zx = DrvROMS(pm, z, 'x'); % XXX-Is this the right type of derivative to take?
    Zy = DrvROMS(pn, z, 'y');
    if (i>1)
        zwt = squeeze(zw(:,:,:,i) - zw(:,:,:,i-1))./ts; % XXX- staggered back 1 timestep...
        zwt = (zwt(:,:,2:end) + zwt(:,:,1:end-1)).*0.5;
%         zwt = 0;
    else
        zwt = 0;
    end
    
    % 3) Load Zonal Momentum Stuff
    U = GetVarROMS(path, 0, {'u', '(1)'}, sliceT, 1);
    Uu = GetVarROMS(path, 0, {'u', '(1)'}, sliceT, 0);
    % Advective Terms
%     Ux = DrvS(pm,z, U, 'x');
%     Uy = DrvS(pn,z, U, 'y');
    Uz = DrvROMS(z, U, 'z');
    
    % XXX-ALTERNATE DERIVS
    Ux = DrvS(pm,z, Uu, 'x', 2);
    Uy = DrvS(pn,z, Uu, 'y', 2);

    % Temporal Difference
    DuDt(:,:,:,i) = (U - UPrev)./ts - Uz.*zwt;
    UPrev = U;
    
    % 4) Load Meridional Velocity
    V = GetVarROMS(path, 0, {'v', '(1)'}, sliceT);
    Vv = GetVarROMS(path, 0, {'v', '(1)'}, sliceT,0);

    % Advective V
%     Vx = DrvS(pm,z, V, 'x');
%     Vy = DrvS(pn,z, V, 'y');
    Vz = DrvROMS(z, V, 'z');
        % XXX-ALTERNATE DERIVS
    % XXX-ALTERNATE DERIVS
   Vx = DrvS(pm,z, Vv, 'x', 3);
    Vy = DrvS(pn,z, Vv, 'y', 3);

     DvDt(:,:,:,i) = (V - VPrev)./ts - Vz.*zwt;
    VPrev = V;
%     Wz = -Ux - Vy;
    
    % 5) Load Buoyancy Terms
    T = GetVarROMS(path, 0, {'temp', '(1)'}, sliceT);
    Tf(:,:,:,i) = T;
    S = GetVarROMS(path, 0, {'salt', '(1)'}, sliceT);
    Sf(:,:,:,i) = S;
%     rho = roms_eos(T, S, 0); % XXX - referencing to 0 db.
%     rho = rho + rho0;
    rho = rho_eos(T, S, z); % CROCO function (checked for consistency 1/12/17)
    Rho(:,:,:,i) = rho;
    B = -g*rho./rho0; % In-Situ B
    % This calculates adiabatically referenced N^2    
    ztemp = permute(z, [3 1 2]);
    zwtemp = permute(squeeze(zw(:,:,1:end-1,i)), [3 1 2]);
    Ttemp = permute(T, [3 1 2]);
    Stemp = permute(S, [3 1 2]);
    [rhoa, bvf] = rho_eos(Ttemp,Stemp,ztemp, zwtemp, g, rho0);
    bvf = permute(bvf, [2 3 1]); % N^2
    

    Bx(:,:,:,i) = DrvS(pm, z, B, 'x');
    By(:,:,:,i) = DrvS(pn, z, B, 'y');
    % - XXX Should be adiabatically leveled first....
%     Bz(:,:,:,i) = DrvROMS(z, B, 'z');
    Bz(:,:,:,i) = bvf;
    % Iterate for vertical integration terms
        DbDt(:,:,:,i) = (B - BPrev)./ts - bvf.*zwt;
    BPrev = B;
    
    W = zwt + U.*Zx + V.*Zy + omega;% see: https://www.myroms.org/forum/viewtopic.php?f=19&t=2139
    
    % Calculate hydro pressure, probably a much better way to do this...
    P = NaN(nx, ny, nz);
    for x =1:nx
        for y =1:ny
             zt = squeeze(z(x,y,:));
            bt = squeeze(B(x,y,:));
            mask = isfinite(bt);
            P(x, y, mask) = -rho0.*B(x, y, end).*Eta(x,y) + rho0.*flipud(cumtrapz(flipud(zt(mask)), flipud(bt(mask))));
        end
    end
    
    % Pressure gradients.
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
% XXX-z shifting  included above
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
xla = 10:nx-10; yla = 10:ny-10; zla = 1:10;
mask = ones(size(Q));
mask = 0*mask;
mask(xla, yla,:,:) = 1;
maskA = mask;
% mask =  (Rho < 1026);
mask(:,:,end,:) = 0;
mask(:,:,1,:) = 0;
mask(1:2,:,:,:) = 0; mask(end-1:end,:,:,:) = 0;
mask(:, 1:2,:,:) = 0; mask(:,end-1:end,:,:) = 0;
dz = diff(zw, 1, 3);
gridvol = abs(repmat(repmat(1./pm, [1 1 nz]).*repmat(1./pn, [1 1 nz]), [1 1 1 nt])).*abs(dz);
vol = squeeze(sum(sum(sum(mask.*gridvol))));

IntegrateQTerms;

%Vertical Terms
dx = repmat(1./pm, [1 1 nt]);
dy = repmat(1./pn, [1 1 nt]);
zlow = zla(1)+1; zhigh = zla(end) - 1;
[JFz_t, dJFz_t    ] = areaIntegrateJVecsROMS(squeeze(JFz(:,:,zhigh,:)), squeeze(mask(:,:,zhigh,:)), dx, dy, ts, vol);
[JFz_b, dJFz_b    ] = areaIntegrateJVecsROMS(squeeze(JFz(:,:,zlow,:)), squeeze(mask(:,:,zlow,:)), dx, dy, ts, vol);
JFzA = JFz_t - JFz_b;
dJFzA = dJFz_t - dJFz_b;
[JBz_t, ~    ] = areaIntegrateJVecsROMS(squeeze(JBz(:,:,zhigh,:)), squeeze(mask(:,:,zhigh,:)), dx, dy, ts, vol);
[JBz_b, ~    ] = areaIntegrateJVecsROMS(squeeze(JBz(:,:,zlow,:)), squeeze(mask(:,:,zlow,:)), dx, dy, ts, vol);
JBzA = JBz_t - JBz_b;
[JAzA, dJAzA] = areaIntegrateJVecsROMS(squeeze(JAz(:,:,zlow,:)), squeeze(maskA(:,:,zlow,:)), dx, dy, ts, vol);
JAzA = -JAzA;
dJAzA = -dJAzA;

%Zonal Terms
xleft = xla(1) ; xright = xla(end) ;
dyl = repmat(squeeze(1./pn(xleft,:).'), [1 nz nt]);
dyr = repmat(squeeze(1./pn(xright,:).'), [1 nz nt]);
% dz = diff(zw,1,3);
% dzl = repmat(squeeze(dz(3,:,:)), [1 1 nt]);
% dzr = repmat(squeeze(dz(end-2,:,:)), [1 1 nt]);
dzl = squeeze(dz(xleft,:,:,:));
dzr = squeeze(dz(xright,:,:,:));
[JFx_l, dJFx_l ] = areaIntegrateJVecsROMS(squeeze(JFx(xleft,:,:,:)), squeeze(mask(xleft,:,:,:)),  dyl, dzl, ts, vol);
[JFx_r, dJFx_r    ] = areaIntegrateJVecsROMS(squeeze(JFx(xright,:,:,:)), squeeze(mask(xright,:,:,:)), dyr, dzr, ts, vol);
JFxA = JFx_r - JFx_l;
dJFxA = dJFx_r - dJFx_l;
[JBx_l, dJBx_l    ] = areaIntegrateJVecsROMS(squeeze(JBx(xleft,:,:,:)), squeeze(mask(xleft,:,:,:)),  dyl, dzl, ts, vol);
[JBx_r, dJBx_r    ] = areaIntegrateJVecsROMS(squeeze(JBx(xright,:,:,:)), squeeze(mask(xright,:,:,:)), dyr, dzr, ts, vol);
JBxA = JBx_r - JBx_l;
dJBxA = dJBx_r - dJBx_l;
[JAx_l, dJAx_l    ] = areaIntegrateJVecsROMS(squeeze(JAx(xleft-1,:,:,:)), squeeze(maskA(xleft-1,:,:,:)),  dyl, dzl, ts, vol);
[JAx_r, dJAx_r    ] = areaIntegrateJVecsROMS(squeeze(JAx(xright+1,:,:,:)), squeeze(maskA(xright+1,:,:,:)), dyr, dzr, ts, vol);
JAxA = JAx_r - JAx_l;
dJAxA = dJAx_r - dJAx_l;
%Meridional Terms
yfront = yla(1); yback = yla(end);
dxf = repmat(squeeze(1./pm(:,yfront)), [1 nz nt]);
dxb = repmat(squeeze(1./pm(:,yback)), [1 nz nt]);
% dzf = repmat(squeeze(dz(:,3,:)), [1 1 nt]);
% dzb = repmat(squeeze(dz(:,end-2,:)), [1 1 nt]);
dzf = squeeze(dz(:,yfront,:,:));
dzb = squeeze(dz(:,yback,:,:));
[JFy_f, ~    ] = areaIntegrateJVecsROMS(squeeze(JFy(:,yfront,:,:)), squeeze(mask(:,yfront,:,:)),  dxf, dzf, ts, vol);
[JFy_b, ~    ] = areaIntegrateJVecsROMS(squeeze(JFy(:,yback,:,:)), squeeze(mask(:,yback,:,:)), dxb, dzb, ts, vol);
JFyA = JFy_b - JFy_f;
[JBy_f, ~    ] = areaIntegrateJVecsROMS(squeeze(JBy(:,yfront,:,:)), squeeze(mask(:,yfront,:,:)),  dxf, dzf, ts, vol);
[JBy_b, ~    ] = areaIntegrateJVecsROMS(squeeze(JBy(:,yback,:,:)), squeeze(mask(:,yback,:,:)), dxb, dzb, ts, vol);
JByA = JBy_b - JBy_f;
[JAy_f, dJAy_f   ] = areaIntegrateJVecsROMS(squeeze(JAy(:,yfront-1,:,:)), squeeze(maskA(:,yfront-1,:,:)),  dxf, dzf, ts, vol);
[JAy_b, dJAy_b    ] = areaIntegrateJVecsROMS(squeeze(JAy(:,yback+1,:,:)), squeeze(maskA(:,yback+1,:,:)), dxb, dzb, ts, vol);
JAyA = JAy_b - JAy_f;
dJAyA = dJAy_b - dJAy_f;
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
% plot(-JFzA-JBzA);
% plot(-(JFzA + JBzA+JAxA+JAyA+JAzA))
% plot(Qas)
hold off
legend('\Delta Q','-J^z_F', '-J_B^z', '-J_F^h', '-J_B^h', '-J_A', 'Sum', 'Location', 'SouthWest')
grid on
xlabel('Timestep (12 hours)');
set(gca, 'FontSize', 18);

%% Q - DIFF

% dQ = Qa + (JFzA + JBzA+JFxA+JFyA+JBxA+JByA+JAxA+JAyA+JAzA);
% % plot(gradient(smooth(dQ.*vol, 20), ts))
% plot(smooth(dQ,20))
% hold on
% plot(-JAxA);
% plot(-JAyA);
% plot(-JAzA);
% plot((JAxA + JAyA));
% hold off
%% %NON CONSERVATIVE Q PLOT

figure
plot(Qa+(JAxA+JAyA+JAzA), 'LineWidth', 2);
hold on
plot(-JFzA);
plot(-JBzA);
plot(-JFxA - JFyA);
plot(-JBxA - JByA);
plot(-(JFzA + JBzA+JFxA+JFyA+JBxA+JByA), '--', 'LineWidth', 2);
% plot(-(JFzA + JBzA+(JFxA+JFyA)+JBxA+JByA), '--', 'LineWidth', 2);

hold off
legend('\Delta Q','-JFz', '-JBz', '-JFh', '-JBh',  'Sum', 'Location', 'SouthWest')
grid on
xlabel('Timestep (12 hours)');
set(gca, 'FontSize', 18);

%%
forcepath = [pardir 'nesea_frc.nc'];
Qor = ncread(forcepath,'shflux'); % Monthly Values W/m^2
EPr = ncread(forcepath,'swflux')./(100.*86400); %Freshwater flux m/s
tau_u = ncread(forcepath,'sustr'); %N/m^2;
tau_v = ncread(forcepath, 'svstr'); %N/m^2
Qor = Qor(xl, yl,:);
EPr = EPr(xl, yl,:);
tau_u = tau_u(xl, yl,:);
tau_v = tau_v(xl, yl,:);

climtime = datenum(2012,1:12, 1);
modeltime = datenum(2012,1,1, 1:12:(nt*12),0,0);
% disp('here')
Qo = NaN(nx, ny, nt);
EP = Qo;
tx = Qo;
ty = Qo;
for x = 1:nx;
    for y = 1:ny
Qo(x, y,:) = interp1(climtime, squeeze(Qor(x,y,:)), modeltime, 'pchip');
EP(x,y,:) = interp1(climtime, squeeze(EPr(x,y,:)), modeltime, 'pchip');
tx(x,y,:) = interp1(climtime, squeeze(tau_u(x,y,:)), modeltime, 'pchip');
ty(x,y,:) = interp1(climtime, squeeze(tau_v(x,y,:)), modeltime, 'pchip');
    end
end

Tmean = squeeze(nanmean(nanmean(nanmean(Tf(:,:,end,:), 4))));
Smean = squeeze(nanmean(nanmean(nanmean(Sf(:,:,end,:), 4))));
%%
bseasonalh =squeeze( nanmean(nanmean(-g*alpha*Qor./(rho0*Cp))));
bseasonals = squeeze(nanmean(nanmean(g*beta*EPr.*35.7)));
bseasonal = bseasonalh+bseasonals;

figure
subplot(2,1,1)
plotyy(1:12, squeeze(nanmean(nanmean(Qor))),1:12, squeeze(nanmean(nanmean(EPr))));
legend('Q_o', 'E-P');
grid on
subplot(2,1,2)
plot(bseasonal, 'LineWidth', 2);
hold on
% plot(bseasonalh);
% plot(bseasonals);
hold off
ylabel('B_o');
grid on
%% GEN THEORY SCALINGS
del = .1;
alpha = 1/rho0.*(rho_eos(Tmean+del, Smean, 0)-rho_eos(Tmean-del, Smean, 0))./(2*del);
beta = 1/rho0.*(rho_eos(Tmean, Smean+del, 0)-rho_eos(Tmean, Smean-del, 0))./(2*del);

Cp = 3994; % XXX - Get ROMS Value
Bo = -g*alpha*Qo./(rho0*Cp) + g*beta*EP.*squeeze(Sf(:,:,end,:));

M2 = squeeze(Bx(:,:,end-2,:).^2 + By(:,:,end-2,:).^2);
H = hkpp; H(:,:,1) = H(:,:,2);
H(H<1) = 1;
% H = repmat(nanmean(nanmean(H)), [nx ny 1]);
JBT = -repmat(f, [1 1 nt]).*Bo./H + 0.16*M2.*H;
JFT = -0.2*M2.*H;

tmag = abs(tx+1i*ty);
de = sqrt(2*.1*H.*sqrt(tmag./rho0)./repmat(f,[1 1 nt]));

JFW = tx./(rho0*H).*squeeze(By(:,:,end-2,:)) - ty./(rho0*H).*squeeze(Bx(:,:,end-2,:));

txa = repmat(nanmean(nanmean(tx)), [nx ny 1]);
tya = repmat(nanmean(nanmean(ty)), [nx ny 1]);

% JFW = txa./(rho0*H).*squeeze(By(:,:,end-2,:)) - tya./(rho0*H).*squeeze(Bx(:,:,end-2,:));
dx = repmat(1./pm, [1 1 nt]);
dy = repmat(1./pn, [1 1 nt]);
[JBTA, dJBTA    ] = areaIntegrateJVecsROMS(JBT, squeeze(mask(:,:,end-1,:)), dx, dy, ts, vol);
[JFTA, dJFTA    ] = areaIntegrateJVecsROMS(JFT, squeeze(mask(:,:,end-1,:)), dx, dy, ts, vol);
[JFWA, dJFWA    ] = areaIntegrateJVecsROMS(JFW, squeeze(mask(:,:,end-1,:)), dx, dy, ts, vol);

%%
figure
plot(Qa+(JAxA+JAyA+JAzA), 'LineWidth', 2);
hold on
plot(-JFzA);
plot(-JBzA);
plot(-JFTA);
plot(-JBTA);
plot(-JFWA)

plot(-(JFTA+JBTA), 'LineWidth', 2, 'LineStyle', '--')
% plot(-(JFzA + JBzA+JFxA + JFyA), 'k','LineWidth', 2)
hold off
legend('\Delta Q_{NC}','-JFz', '-JBz', '-JF_{SCALE}', '-JB_{SCALE}', '-JW_{SCALE}', 'SUM-SCALING','Location', 'SouthWest')
grid on

%% SCATTER PLOT
subplot(2,1,1)
scatter(-JFzA, -(JFTA))
cr = corr(JFzA, JFTA);
title(num2str(cr));
grid on
subplot(2,1,2)
scatter(-JBzA, -JBTA)
grid on
%% NON-CONS d/dt terms.
s = 3;
figure
plot(Qt);
hold on
plot(-dJBTA);
plot(-dJFTA);
plot(-(dJBTA+dJFTA));
% plot(-smooth(dJBTA+dJFTA,s));

hold off
grid on
%%

plot(-JAzA);
hold on
plot(-JAxA);
plot(-JAyA);
plot(-JAxA-JAyA-JAzA)
hold off