% NOTING all things to be looked at later:
% XXX - grid is spherical.
% XXX - Probably should directly calculate Pressure gradients (not press)

%PARAMETERS TO CHANGE'
% pardir = '/groups/thomas1/jacob13/GulfStream/NESEB/';
% basepath = [pardir 'HIS/'];
% path1 = '/groups/thomas1/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc';
% ntp = 2;
% nt = 400;
% inidate = datenum(2012, 1, 1); %XXX- assuming it starts on Jan 1...?
% dateoff = datenum(0, 0, 40);
% modeltime = inidate + dateoff + datenum(0, 0, 0, 1:1:nt, 0,0);
% forcepath = [pardir 'nesea_frc.nc'];

pardir = '/data/thomas/jacob13/GULFZ/';
basepath = [pardir 'HIS/'];
path1 = [pardir 'gulfz_grd.nc'];
forcepath = [pardir 'gulfz_frc.nc'];
ntp = 5;
nt = 100;345;
offset = 160;
modeltime = datenum(2012,8,(26+offset):1:(26+offset+(nt)-1), 0,0,0);

%
files = dir([basepath,'*.nc']);
zl = 1:50;
nz = length(zl);
xl = 700:1000;800:850;
yl = 050:600;400:450;350:400;
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

latst = ncread(path1, 'lat_rho');
lats = latst(xl, yl);
lonst = ncread(path1, 'lon_rho');
lons = lonst(xl, yl);
disp(['Sub-domain borders (clockwise from bottom left): (', num2str(lats(1,1)), ',',num2str(lons(1,1)),')  (',...
    num2str(lats(1,end)),',',num2str(lons(1,end)),')  (', num2str(lats(end,end)),',',num2str(lons(end,end)),')  (',...
    num2str(lats(end,1)),',',num2str(lons(end,1)),')']);
% FULL DOMAIN
path = [basepath, files(floor(offset./ntp)).name];

T = GetVarROMS(path, 0, {'temp', '(1)'}, {0, 0,0, 0});

subplot(2,1,1)
% pcolor( squeeze(lonst(:,1)),squeeze(latst(1000,:)),(squeeze(T(:,:,end,1))).'); shading interp
pcolor( (squeeze(T(:,:,end,1))).'); shading interp

hold on
rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
hold off
colorbar;


%%
% pnmetric = repmat(pn, [1 1 nz nt]);
% Ts = NaN(1602, 922, nt); 
% Hs = Ts;
% for i=1:50;length(files);
%     
%     disp([num2str(i), '/', num2str((length(files)-2)*2)]);
%     fileind = ceil((i+offset)/ntp);
%     name = files(fileind).name;
%     path = [basepath, name];
%     sliceind = mod(i,ntp);
%     if sliceind==0; sliceind=5;end
%     disp(['File Ind: ', num2str(fileind), '  Slice Ind: ',num2str(sliceind)])
%     sliceT = {0, 0, [50 50],[sliceind sliceind]};
%     Ts(:,:,i) = squeeze(GetVarROMS(path, 0, {'temp', '(1)'}, sliceT));
%     Hs(:,:,i) = squeeze(GetVarROMS(path, 0, {'hbls', '(1)'}, sliceT));
% end
% 
% %%
% figure
% for i=1:nt
%    pcolor(squeeze(Ts(:,:,i)).'); shading interp
%    colorbar
%    hold on
% rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
% hold off
% title(num2str(i))
% drawnow
%  pause(0.001);
% end
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
% offset = 0; % Sets the 'start date'
[nx ny] = size(pm);
UPrev = NaN;
VPrev = NaN;
BPrev = NaN;
DuDt = NaN(nx, ny, nz, nt);

Q = DuDt;
Bx = DuDt;
By = DuDt;
Bz = DuDt;
Tf = DuDt;
Uf = DuDt;
Vf = DuDt;
Wf = DuDt;
Sf = DuDt;
Rho = DuDt;
zw = NaN(nx, ny, nz+1, nt);
hkpp = NaN(nx, ny, nt);
OMEGAX = DuDt;
OMEGAY = DuDt;
OMEGAZ = DuDt;
JAz = DuDt; JAx = DuDt; JAy = DuDt;
dADV = JAz;

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
    Uz = DrvROMS(z, U, 'z');
    Uy = DrvS(pn,z, Uu, 'y', 2);
    
    
    % 4) Load Meridional Velocity
    V = GetVarROMS(path, 0, {'v', '(1)'}, sliceT);
    Vv = GetVarROMS(path, 0, {'v', '(1)'}, sliceT,0);

    % Advective V

    Vz = DrvROMS(z, V, 'z');
   Vx = DrvS(pm,z, Vv, 'x', 3);

    
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
    Bz(:,:,:,i) = DrvROMS(z, B, 'z');
%     Bz(:,:,:,i) = bvf;
    % Iterate for vertical integration terms
    
    W = zwt + U.*Zx + V.*Zy + omega;% see: https://www.myroms.org/forum/viewtopic.php?f=19&t=2139
    Wy = DrvS(pn, z, W, 'y');
    Wx = DrvS(pm, z, W, 'x');
    
    Uf(:,:,:,i) = U; Vf(:,:,:,i) = V; Wf(:,:,:,i) = W;
    OMEGAX(:,:,:,i) = - Vz;
    OMEGAY(:,:,:,i) = + Uz;
    OMEGAZ(:,:,:,i) = repmat(f, [1 1 nz]) + Vx - Uy;
    
    Q(:,:,:,i) = OMEGAX(:,:,:,i).*Bx(:,:,:,i) + OMEGAY(:,:,:,i).*By(:,:,:,i)+OMEGAZ(:,:,:,i).*Bz(:,:,:,i);
    
  % SHOULD KEEP ADVECTIVE TERMS...  
      JAz(:,:,:,i) = W.*Q(:,:,:,i); JAx(:,:,:,i) = U.*Q(:,:,:,i); JAy(:,:,:,i) = V.*Q(:,:,:,i);
      
%       JAx(2:end,:,:,i) = Uu(1:end,:,:).*(Q(1:end-1,:,:,i) + Q(2:end,:,:,i))./2;
%       JAx(2:end-1,:,:,i) = (JAx(2:end-1,:,:,i) + JAx(3:end,:,:,i))./2; 
%       
%         JAy(:,2:end,:,i) = Vv(:,1:end,:).*(Q(:,1:end-1,:,i) + Q(:,2:end,:,i))./2;
%       JAy(:,2:end-1,:,i) = (JAy(:,2:end-1,:,i) + JAy(:,3:end,:,i))./2; 
%     dADV(:,:,:,i) = DrvROMS(z, JAz(:,:,:,i), 'z') + DrvS(pm,z, JAx(:,:,:,i), 'x') + DrvS(pn,z,JAy(:,:,:,i),'y');
    dADV(:,:,:,i) = W.*DrvROMS(z, Q(:,:,:,i), 'z') + U.*DrvS(pm,z, Q(:,:,:,i), 'x') + V.*DrvS(pn,z,Q(:,:,:,i),'y');
end


%%
xla = 2:nx-1; yla = 2:ny-1; zla =1:50;
mask = zeros(size(Q));
% mask = 0*mask;
mask(:, :,zla,:) = 1;
% mask = (Tf>18) & (Tf<20); % Mode Water Mask
% mask = (Tf>18) & (Tf<20);
% mask = (Rho> 1025.6) & (Rho<1026.2);
mask = (Rho>1025.8) & (Rho<1026.2);
% mask = (rho>1027)& (Rho<1028);
mh = hkpp>100;
mh = permute(repmat(mh, [1 1 1 nz]), [1 2 4 3]);
% mask = mask.*mh;
maskA = mask;
% mask =  (Rho < 1026);
% mask(:,:,end,:) = 0;
% mask(:,:,1,:) = 0;
% mask(1:2,:,:,:) = 0; mask(end-1:end,:,:,:) = 0;
% mask(:, 1:2,:,:) = 0; mask(:,end-1:end,:,:) = 0;
dz = diff(zw, 1, 3);
gridvol = abs(repmat(repmat(1./pm, [1 1 nz]).*repmat(1./pn, [1 1 nz]), [1 1 1 nt])).*abs(dz);
vol = squeeze(sum(sum(sum(mask.*gridvol))));
vol = nanmean(vol);

IntegrateQTerms
disp('here')
%Vertical Terms
dx = repmat(1./pm, [1 1 nt]);
dy = repmat(1./pn, [1 1 nt]);
zlow = zla(1); zhigh = zla(end) ;
% JAzt = NaN(size(JAz)); 
% JAzt(:,:,1:end-1,:) = (JAz(:,:,1:end-1,:) + JAz(:,:,2:end,:))./2;
JAzt = JAz;
[JAzAb, dJAzAb] = areaIntegrateJVecsROMS(squeeze(JAzt(:,:,zlow,:)), squeeze(maskA(:,:,zlow,:)), dx, dy, ts, vol);
[JAzAt, dJAzAt] = areaIntegrateJVecsROMS(squeeze(JAzt(:,:,zhigh,:)), squeeze(maskA(:,:,zhigh,:)), dx, dy, ts, vol);

JAzA = JAzAt-JAzAb;
dJAzA = dJAzAt-dJAzAb;

%Zonal Terms
xleft = xla(1) ; xright = xla(end) ;
dyl = repmat(squeeze(1./pn(xleft,:).'), [1 nz nt]);
dyr = repmat(squeeze(1./pn(xright,:).'), [1 nz nt]);
dzl = squeeze(dz(xleft,:,:,:));
dzr = squeeze(dz(xright,:,:,:));
[JAx_l, dJAx_l    ] = areaIntegrateJVecsROMS(squeeze(JAx(xleft,:,:,:)), squeeze(maskA(xleft,:,:,:)),  dyl, dzl, ts, vol);
[JAx_r, dJAx_r    ] = areaIntegrateJVecsROMS(squeeze(JAx(xright,:,:,:)), squeeze(maskA(xright,:,:,:)), dyr, dzr, ts, vol);
JAxA = JAx_r - JAx_l;
dJAxA = dJAx_r - dJAx_l;

yfront = yla(1); yback = yla(end);
dxf = repmat(squeeze(1./pm(:,yfront)), [1 nz nt]);
dxb = repmat(squeeze(1./pm(:,yback)), [1 nz nt]);
dzf = squeeze(dz(:,yfront,:,:));
dzb = squeeze(dz(:,yback,:,:));
[JAy_f, dJAy_f   ] = areaIntegrateJVecsROMS(squeeze(JAy(:,yfront,:,:)), squeeze(maskA(:,yfront,:,:)),  dxf, dzf, ts, vol);
[JAy_b, dJAy_b    ] = areaIntegrateJVecsROMS(squeeze(JAy(:,yback,:,:)), squeeze(maskA(:,yback,:,:)), dxb, dzb, ts, vol);
JAyA = JAy_b - JAy_f;
dJAyA = dJAy_b - dJAy_f;

QADV = JAyA + JAxA + JAzA;
QADV = QADV - QADV(1);
dQADV = dJAyA + dJAxA + dJAzA;
%%

figure
subplot(2,1,1)
plot(Qa, 'LineWidth', 2);
hold on
plot(-QADV);
plot(Qa + QADV, '--');
% plot(-JFzA);
% plot(-JBzA);
% plot(-JFxA - JFyA);
% plot(-JBxA - JByA);
% plot(-JAxA - JAyA-JAzA);
% plot(-(JFzA + JBzA+JFxA+JFyA+JBxA+JByA+JAxA+JAyA+JAzA), '--', 'LineWidth', 2);
% plot(-JFzA-JBzA);
% plot(-(JFzA + JBzA+JAxA+JAyA+JAzA))
% plot(Qas)
hold off
legend('\Delta Q','-J^z_F', '-J_B^z', '-J_F^h', '-J_B^h', '-J_A', 'Sum', 'Location', 'SouthWest')
grid on
xlabel('Timestep (12 hours)');
set(gca, 'FontSize', 18);

subplot(2,1,2)
plot(squeeze(nanmean(nanmean(mask(:,:,end-1,:)))))
hold on
plot(squeeze(nanmean(nanmean(mask(:,:,2,:)))), '--');
hold off

%%

figure
subplot(2,1,1)
plot(Qt, 'LineWidth', 2);
hold on
plot(-dQADV);
plot(Qt + dQADV, '--');
% plot(-JFzA);
% plot(-JBzA);
% plot(-JFxA - JFyA);
% plot(-JBxA - JByA);
% plot(-JAxA - JAyA-JAzA);
% plot(-(JFzA + JBzA+JFxA+JFyA+JBxA+JByA+JAxA+JAyA+JAzA), '--', 'LineWidth', 2);
% plot(-JFzA-JBzA);
% plot(-(JFzA + JBzA+JAxA+JAyA+JAzA))
% plot(Qas)
hold off
legend('\Delta Q','-J^z_F', '-J_B^z', '-J_F^h', '-J_B^h', '-J_A', 'Sum', 'Location', 'SouthWest')
grid on
xlabel('Timestep (12 hours)');
set(gca, 'FontSize', 18);

subplot(2,1,2)
plot(squeeze(nanmean(nanmean(mask(:,:,end-1,:)))))
hold on
plot(squeeze(nanmean(nanmean(mask(:,:,2,:)))), '--');
hold off
%%
if true
    zeta = OMEGAZ - repmat(f, [1 1 nz nt]);
    cl = [1024 1027];
    for i=1:nt
        subplot(2,1,1)
       pcolor(squeeze(Rho(:,:,end-1,i)).'); shading interp

       colorbar
       hold on
       contour(squeeze(mask(:,:,end-1,i)).', [0 1], 'k');
       hold off
       set(gca, 'clim', cl);
       subplot(2,1,2)
       pcolor((squeeze(zeta(:,:,end-1,i))./f).'); shading interp
          hold on
       contour(squeeze(mask(:,:,end-1,i)).', [0 1], 'k');
       hold off
       set(gca, 'clim', [-2 2]);
       colorbar
       pause(0.001)
    end
end

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

% figure
% plot(Qa+(JAxA+JAyA+JAzA), 'LineWidth', 2);
% hold on
% plot(-JFzA);
% plot(-JBzA);
% plot(-JFxA - JFyA);
% plot(-JBxA - JByA);
% plot(-(JFzA + JBzA+JFxA+JFyA+JBxA+JByA), '--', 'LineWidth', 2);
% % plot(-(JFzA + JBzA+(JFxA+JFyA)+JBxA+JByA), '--', 'LineWidth', 2);
% 
% hold off
% legend('\Delta Q','-JFz', '-JBz', '-JFh', '-JBh',  'Sum', 'Location', 'SouthWest')
% grid on
% xlabel('Timestep (12 hours)');
% set(gca, 'FontSize', 18);

%%

Qor = ncread(forcepath,'shflux'); % Monthly Values W/m^2
EPr = ncread(forcepath,'swflux')./(100.*86400); %Freshwater flux m/s
tau_u = ncread(forcepath,'sustr'); %N/m^2;
tau_v = ncread(forcepath, 'svstr'); %N/m^2
Qor = Qor(xl, yl,:);
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
climtime = datenum(2012,1, 15:30:(360*2));

% disp('here')
Qo = NaN(nx, ny, nt);
EP = Qo;
tx = Qo;
sst = Qo;
ty = Qo;
for x = 1:nx;
    disp(num2str(x))
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
%%
% mlmask = Rho > repmat(Rho(:,:,end,:), [1 1 nz 1])+0.125;
% mlmask = double(mlmask);
% mlmask(mlmask==0) = NaN;
% zf = diff(zw, 1, 3);
% zf = -permute(repmat(zm, [1 nx ny nt]), [2 3 1 4]);
% mldepth =squeeze( min(zf.*mlmask, [], 3));
%% GEN THEORY SCALINGS
del = .1;

rhot = rho_eos(squeeze(Tf(:,:,end,:)), squeeze(Sf(:,:,end,:)), 0);
alpha = 1./rhot.*(rho_eos(squeeze(Tf(:,:,end,:))+del, squeeze(Sf(:,:,end,:)), 0)-rho_eos(squeeze(Tf(:,:,end,:))-del, squeeze(Sf(:,:,end,:)), 0))./(2*del);
beta = 1./rhot.*(rho_eos(squeeze(Tf(:,:,end,:)), squeeze(Sf(:,:,end,:))+del, 0)-rho_eos(squeeze(Tf(:,:,end,:)), squeeze(Sf(:,:,end,:))-del, 0))./(2*del);

Cp = 3994; % XXX - Get ROMS Value
ssta =  squeeze(Tf(:,:,end,:))-sst;
Qeff = Qo - 30*ssta;
Bo = g*alpha.*(Qeff)./(rho0*Cp) +  g*beta.*EP.*squeeze(Sf(:,:,end,:));
% Bo = g*alpha.*(Qo)./(rho0*Cp) +  g*beta.*EP.*squeeze(Sf(:,:,end,:));
% zf = repmat(zw, [1 1 1 nt]);
zf = zw(:,:,1:end-1,:);
H = hkpp; H(:,:,1) = H(:,:,2);
% H(mldepth>H) = mldepth(mldepth>H);
% H = squeeze(mldepth);
% H(H<12) = 12;
ms = double(zf > -permute(repmat(H, [1 1 1 nz]), [1 2 4 3]));
ms(ms ==0) = NaN;
M2 = squeeze(squeeze(nanmean(ms.*Bx(:,:,:,:), 3)).^2 + squeeze(nanmean(ms.*By(:,:,:,:), 3)).^2); %  XX- not averaging over BL...
m = ~isfinite(M2);  
Ma = squeeze(Bx(:,:,end-1,:)).^2 + squeeze(By(:,:,end-1,:)).^2;
M2(m) = Ma(m); % Use uppermost gradient value if H is too small...


tmag = abs(tx+1i*ty);
de = sqrt(2*.1*H.*sqrt(tmag./rho0)./repmat(f,[1 1 nt]));


H(H<6) = 6;
% H = repmat(nanmean(nanmean(H)), [nx ny 1]);
% JBT = repmat(f, [1 1 nt]).*Bo./H + 0.15*M2.*H;
% JBT(hkpp < 20) = NaN;
% JBTs = JBT - 0.15.*M2.*H;
JFT = -0.2*M2.*H;
% H(H<24) = 24;
Hm =H;
Hm(Hm<23 & Qeff>0) = 23; % Logic here is if the boundary layer is shallow and flux is into ML, correct scale is absorption scale.
JBTs = repmat(f, [1 1 nt]).*Bo./Hm;
% JBTs(JBTs<0) = 0;
JBTe = 0.15.*M2.*H;
JBT = JBTs + JBTe;

disp('here')

R = 0.05.*(H.^2.*M2)./(Bo.*repmat(f, [1 1 nt]));
% H(H<de) = de(H<de);
% H = H;
JFW = -tx./(rho0*H).*squeeze(By(:,:,end-2,:)) + ty./(rho0*H).*squeeze(Bx(:,:,end-2,:));

txa = repmat(nanmean(nanmean(tx)), [nx ny 1]);
tya = repmat(nanmean(nanmean(ty)), [nx ny 1]);
maskh = mask;
% maskh = mask.*mh;
% JFW = txa./(rho0*H).*squeeze(By(:,:,end-2,:)) - tya./(rho0*H).*squeeze(Bx(:,:,end-2,:));
dx = repmat(1./pm, [1 1 nt]);
dy = repmat(1./pn, [1 1 nt]);
[JBTA, dJBTA    ] = areaIntegrateJVecsROMS(JBT, squeeze(maskh(:,:,end-1,:)), dx, dy, ts, vol);
[JFTA, dJFTA    ] = areaIntegrateJVecsROMS(JFT, squeeze(maskh(:,:,end-1,:)), dx, dy, ts, vol);
[JFWA, dJFWA    ] = areaIntegrateJVecsROMS(JFW, squeeze(maskh(:,:,end-1,:)), dx, dy, ts, vol);
[JBTAs, dJBTAs    ] = areaIntegrateJVecsROMS(JBTs, squeeze(maskh(:,:,end-1,:)), dx, dy, ts, vol);
[JBTAe, dJBTAe ] =  areaIntegrateJVecsROMS(JBTe, squeeze(maskh(:,:,end-1,:)), dx, dy, ts, vol);

%%
figure
subplot(3,1,1)
plot(Qa.*vol  , 'LineWidth', 2);
hold on
plot(smooth((Qa +QADV).*vol, 2), 'LineWidth', 2);

plot(-JFTA.*vol);
plot(-JBTA.*vol);
plot(-JFWA.*vol)
plot(-JBTAs.*vol, '--');
plot(-JBTAe.*vol, '--');

plot(-(JFTA+JBTAs+JBTAe+JFWA).*vol, 'LineWidth', 2, 'LineStyle', '--')
% plot(-(JFzA + JBzA+JFxA + JFyA), 'k','LineWidth', 2)
hold off
legend('\Delta Q_{NC}','Smoothed','-JF_{SCALE}', '-JB_{SCALE}', '-JW_{SCALE}','JD_S', 'JD_E', 'SUM-SCALING','Location', 'SouthWest')
grid on

subplot(3,1,2)
plotyy(1:nt, squeeze(nanmean(nanmean(H.*squeeze(mask(:,:,end-1,:))))), 1:nt, squeeze(nanmean(nanmean(Qo.*squeeze(mask(:,:,end-1,:))))));
grid on

subplot(3,1,3)
s=3;
hold on
plot(-smooth(dJBTA, s));
plot(-smooth(dJFTA,s));
plot(-smooth(dJFWA,s));
plot(-smooth(dJBTAs*1+ ((dJBTA -dJBTAs)*1 +dJFTA)+dJFWA,s));
% plot(-smooth(dJBTAs,s))
% plot(-smooth(dJBTA+dJFTA,s));

hold off
grid on
legend('J_D', 'J_F', 'JF_W', 'Sum')
%%
st =150;
qnc = (Qa+QADV).*vol;
plot(qnc(st:end) - qnc(st));
hold on
jt = -(JFTA + JBTA + JFWA).*vol;
plot(jt(st:end) - jt(st));
hold off
%% SCATTER PLOT
subplot(2,1,1)
scatter(-JFzA, -(JFTA))
cr = corr(JFzA, JFTA);
title(num2str(cr));
onetoone
grid on
subplot(2,1,2)
scatter(-JBzA, -JBTA)
onetoone
grid on
%%
scatter(Qa + JFTA+JBTA+JFWA, -QADV);
corr(Qa + JFTA+JBTA+JFWA, -QADV)

%% NON-CONS d/dt terms.
s =12;
Qta = gradient((Qa+QADV).*vol, ts);
Qtana = gradient(Qa*vol, ts);
% Qta = smooth(Qt + dQADV, s);
figure
plot(smooth(Qta,s));
hold on
plot(-smooth(dJBTA, s));
plot(-smooth(dJFTA,s));
plot(-smooth(dJFWA,s));
plot(-smooth(dJBTAs+ ((dJBTA -dJBTAs)*1 +dJFTA)+dJFWA,s));
plot(-smooth(dJBTAs,s))
plot(smooth(Qtana,s))
% plot(-smooth(dJBTA+dJFTA,s));

hold off
grid on
legend('dQdt', 'J_D', 'J_F', 'JF_W', 'Sum')
%%

plot(-JAzA);
hold on
plot(-JAxA);
plot(-JAyA);
plot(-JAxA-JAyA-JAzA)
hold off
legend('z', 'x', 'y', 'sum')
%%

plotyy(1:350, Qa , 1:350, JAxA);

%%
scatter(QADV, squeeze(nanmean(nanmean(R.*squeeze(mask(:,:,end-1,:))))))