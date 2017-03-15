% NOTING all things to be looked at later:
% XXX - grid is spherical.
% XXX - Probably should directly calculate Pressure gradients (not press)
tic

%%
% CAN REPRODUCE THE ISSUE IN NESEA using nt = 43; offset = 150; xl =
% 800:1050; yl = 300:500. Working theory is advection of GS meander into
% domain is not captured by fluxes at edge of isopycnal.

%PARAMETERS TO CHANGE'
pardir = '/groups/thomas1/jacob13/GulfStream/NESEA/';
basepath = [pardir 'HIS/'];
path1 = '/groups/thomas1/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc';
ntp = 2;
nt = 195;
% nt = 400;
% nt = 210;
% nt = 43;
offset = 0;
% offset = 81; % use for NESEB comparison.
% nt = 35; % use for NESB comparison.

% offset = 0; % Use for NESEB
% nt = 418; % Use for NESEB

forcepath = [pardir 'nesea_frc.nc'];
bdom = false;

% % % 
% pardir = '/data/thomas/jacob13/GULFZ/';
% basepath = [pardir 'HIS/'];
% path1 = [pardir 'gulfz_grd.nc'];
% forcepath = [pardir 'gulfz_frc.nc'];
% ntp = 5;
% bdom = true;
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

zl = 1:50;
nz = length(zl);
xl =  1:2002;100:300;1000:1150;1200:1300;
yl =  1:1602; 1:150;250:350;

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
pm = ncread(path, 'pm');
pm = pm(xl, yl);
% pmmetric = repmat(pm, [1 1 nz nt]);

pn = ncread(path, 'pn');
pn = pn(xl,yl);
dx = 1./pm;
dy = 1./pn;
dxf = repmat(dx, [1 1 nz]);
dyf = repmat(dy, [1 1 nz]);

h = ncread(path, 'h');

latst = ncread(path1, 'lat_rho');
lats = latst(xl, yl);
lonst = ncread(path1, 'lon_rho');
lons = lonst(xl, yl);
disp(['Sub-domain borders (clockwise from bottom left): (', num2str(lats(1,1)), ',',num2str(lons(1,1)),')  (',...
    num2str(lats(1,end)),',',num2str(lons(1,end)),')  (', num2str(lats(end,end)),',',num2str(lons(end,end)),')  (',...
    num2str(lats(end,1)),',',num2str(lons(end,1)),')']);
% FULL DOMAIN
path = [basepath, files(ceil((nt+offset)./ntp)).name];
 
%% PLOTS

% subplot(2,1,2)
figure
pcolor(h.'); shading interp
hold on
rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
hold off
colorbar

T = GetVarROMS(path, 0, {'temp', '(1)'}, {0, 0,0, 0});

% subplot(2,1,1)

% pcolor( (squeeze(T(:,:,end,1))).'); shading interp
% pcolor( squeeze(lonst(:,1)),squee`ze(latst(1000,:)),(squeeze(T(:,:,end,1))).'); shading interp
% axesm mercator; framem; gridm; axis on; tightmap
% close all
figure
axesm('mercator', 'MapLatLimit', [26.75 47], 'MapLonLimit', [-90 -55],...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5, 'PLineLocation', 5, 'MLineLocation', 5);
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(latst(:,:)),squeeze(lonst(:,:)),(squeeze(T(:,:,end,1)))); shading interp

lonf = [lons(1,:) lons(:,end).' fliplr(lons(end,:)) flipud(lons(:,1)).'];
latf = [lats(1,:) lats(:,end).' fliplr(lats(end,:)) flipud(lats(:,1)).'];
hold on
hf = fillm(latf, lonf, 'k');
set(hf, 'FaceColor', 'none');
set(hf, 'LineWidth', 2);
% set(h, 
% rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
hold off
% colorbar;


%% LOAD FLUXES

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
% climtime = datenum(2012,1, 15:30:(360*2));
ds = 86400;
climtime = ((360+15)*ds):30*ds:(360*2*ds+360*ds);
% disp('here')
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
% 
% for x = 1:nx;
%     disp(num2str(x))
%     for y = 1:ny
%         try
%         Qo(x, y,:) = interp1(climtime, squeeze(Qor(x,y,:)), modeltime, 'pchip', 0);
%         EP(x,y,:) = interp1(climtime, squeeze(EPr(x,y,:)), modeltime, 'pchip', 0);
%         tx(x,y,:) = interp1(climtime, squeeze(tau_u(x,y,:)), modeltime, 'pchip', 0);
%         ty(x,y,:) = interp1(climtime, squeeze(tau_v(x,y,:)), modeltime, 'pchip', 0);
%         sst(x,y,:) = interp1(climtime, squeeze(SST(x,y,:)), modeltime, 'pchip', 0);
%         SW(x,y,:) = interp1(climtime, squeeze(SWR(x,y,:)), modeltime, 'pchip', 0);
%         catch ME
%             disp(ME)
%         end
%     end
% end


tmag = abs(tx+1i.*ty);

% clear Qor EPr tau_u tau_v SST SWR
toc


%% Calculate Flux Terms
%need du/dt + ADV + CORI + PRESS
% offset = 0; % Sets the 'start date'
[nx ny] = size(pm);

DuDt = single(NaN(nx, ny, nz));
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
zw = single(NaN(nx, ny, nz+1));
hkpp = single(NaN(nx, ny));
OMEGAX = DuDt;
OMEGAY = DuDt;
OMEGAZ = DuDt;
JAz = DuDt; JAx = DuDt; JAy = DuDt;
% dADV = JAz;

name = files(1).name;
path = [basepath, name];
theta_s = ncreadatt(path, '/', 'theta_s');
theta_b = ncreadatt(path, '/', 'theta_b');
hc = ncreadatt(path, '/', 'hc');
% if (size(h)> size([xl yl]))
h = h(xl, yl);
% end
mask = zeros(nx, ny, nz, nt);
vol = NaN(nt, 1);
Qa = NaN(nt, 1);
dJAzA = Qa;
dJAxA = Qa;
dJAyA = Qa;
Qat = Qa;
dJBTA = Qa;
dJFTA = Qa;
dJFWA = Qa;
dJBTAs = Qa;
dJBTAe = Qa;
dJBFWe = Qa;
dJENTa = Qa;
hm = Qa;
Qm = hm;
% Bo = hm;
omegazs = NaN(nx, ny, nt);
% Qf = NaN(nx, ny, nz, nt);

parfor i=1:nt;
    % 1) Housekeeping
    disp([num2str(i), '/', num2str((length(files)-ntp)*ntp)]);
    fileind = ceil((i+offset)/ntp);
    name = files(fileind).name;
    path = [basepath, name];
    sliceind = mod(i,ntp);
    if sliceind==0; sliceind=ntp;end
    sliceT = {slice{1}, slice{2}, slice{3},[sliceind sliceind]};


        outstruct = ROMSFASTPV(path, sliceT, pm, pn, squeeze(sst(:,:,i)), squeeze(Qo(:,:,i)), squeeze(EP(:,:,i)),...
            squeeze(SW(:,:,i)), rho0, Cp, squeeze(tx(:,:,i)), squeeze(ty(:,:,i)), squeeze(tmag(:,:,i)), dxf, dyf,...
            theta_s, theta_b, hc, h, zl, g, f, zmin);
        
        
        mask(:,:,:,i) = outstruct.mask;
        vol(i) = outstruct.vol;
        Qat(i) = outstruct.Qat;
        dJAzA(i) = outstruct.dJAzA;
        dJAxA(i) = outstruct.dJAxA;
        dJAyA(i) = outstruct.dJAyA;
        dJBTA(i) = outstruct.dJBTA;
        dJFTA(i) = outstruct.dJFTA;
        dJFWA(i) = outstruct.dJFWA;
        dJBTAs(i) = outstruct.dJBTAs;
        dJBTAe(i) = outstruct.dJBTAe;
        dJBFWe(i) = outstruct.dJBFWe;
        dJENTa(i) = outstruct.dJENTa;
        hm(i) = outstruct.hm;
        Qm(i) = outstruct.Qm;
%         Bo(i) = outstruct.Bo;
        omegazs(:,:,i) = outstruct.omegazs;
%         Qf(:,:,:,i) = outstruct.Q;
%         Tf(:,:,:,i) = outstruct.T;
end
        disp(['Elapsed Time: ', num2str(toc./60)]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do post calculations
% vol = nanmean(vol);
Qt = gradient(Qat, ts); %Take the time derivative

% Try to work around volume issues.
m = smooth(abs(gradient(vol)./vol),1) < 0.05;
Qa = NaN(length(Qt),1);
Qa(m) = cumtrapz(modeltime(m).',Qt(m));
Qa = interp1(modeltime(m), Qa(m), modeltime).'./vol;

% 'Correct' way to do it.
Qa = Qat-Qat(1); % Interested in Delta Q
Qa = Qa./vol; % Average PV substance in the volume (just a normalization factor).


QADV = cumtrapz(dJAzA + dJAxA + dJAyA).*ts./vol;

% QADV(m) = cumtrapz(modeltime(m), dJAzA(m) + dJAxA(m) + dJAyA(m));
% QADV = interp1(modeltime(m), QADV(m), modeltime).'./vol;

JBTA = cumtrapz(dJBTA).*ts./vol;
JFTA = cumtrapz(dJFTA).*ts./vol;
JFWA = cumtrapz(dJFWA).*ts./vol;
JBTAs = cumtrapz(dJBTAs).*ts./vol;
JBTAe = cumtrapz(dJBTAe).*ts./vol;
JFWAe = cumtrapz(dJBFWe).*ts./vol;
JENTA = cumtrapz(dJENTa).*ts./vol;

toc

%%
figure
% subplot(3,1,1)
plot(modeltime, (Qa+QADV).*vol  , 'LineWidth', 2);
hold on
% plot(smooth((Qa +QADV).*vol, 2), 'LineWidth', 2);

plot(modeltime, -JFTA.*vol, 'LineWidth', 2);
plot(modeltime, -JBTA.*vol, 'LineWidth', 2);
plot(modeltime, -JFWA.*vol, 'LineWidth', 2)
% plot(modeltime, -JBTAs.*vol, '--');
% plot(modeltime, -JBTAe.*vol, '--');


plot(modeltime, -(JFTA+JBTAs+JBTAe+JFWA+JENTA).*vol, 'LineWidth', 2, 'LineStyle', '--')
plot(modeltime, -JENTA.*vol, '--');
% plot(-(JFzA + JBzA+JFxA + JFyA), 'k','LineWidth', 2)
hold off
legend('\Delta \int_V Q + \int_t\int_\Sigma uQ \cdot n','-\int_t\int_A J_{F_{GEO}}', '-\int_t\int_A J_{D}', '-\int_t\int_A J_{F_{WIND}}', '-\int_t\int_A (J_{F_{GEO}} + J_D + J_{F_{WIND}})','Location', 'NorthWest')
grid on
title(['1025.9 < \rho < 1026.3'])

ylabel('m^3 s^{-3}');
set(gca, 'FontSize', 16);
datetick('x')
set(gca, 'xlim', [modeltime(1) modeltime(end)]);
set(gcf, 'Color', 'w');
% subplot(3,1,2)
% plotyy(1:nt, squeeze(nanmean(nanmean(H.*squeeze(mask(:,:,end-1,:))))), 1:nt, squeeze(nanmean(nanmean(Qo.*squeeze(mask(:,:,end-1,:))))));
% grid on
% 
% subplot(3,1,3)
% s=3;
% hold on
% plot(-smooth(dJBTA, s));
% plot(-smooth(dJFTA,s));
% plot(-smooth(dJFWA,s));
% plot(-smooth(dJBTAs*1+ ((dJBTA -dJBTAs)*1 +dJFTA)+dJFWA,s));
% % plot(-smooth(dJBTAs,s))
% % plot(-smooth(dJBTA+dJFTA,s));
% 
% hold off
% grid on
% legend('J_D', 'J_F', 'JF_W', 'Sum')
%%

%% SCATTER PLOT

%%
%% NON-CONS d/dt terms.
s =1;
Qta = gradient((Qa+QADV).*vol, ts);

Qtana = gradient(Qa.*vol, ts);
% Qta = smooth(Qt + dQADV, s);
figure
subplot(2,1,1)
plot(smooth(Qta,s), 'LineWidth', 2);
hold on
plot(-smooth(dJBTAe  +dJBTAs +dJFTA+dJFWA+dJENTa,s), 'LineWidth', 2);
plot(-smooth(dJBTA, s));
plot(-smooth(dJFTA,s));
plot(-smooth(dJFWA,s));

plot(-smooth(dJBTAs,s),'--')
plot(-smooth(dJENTa,s))
% plot(-smooth(dJBTA+dJFTA,s));
hold off
grid on
legend('\partial/\partialt \int_V Q + \int_\Sigma uQ \cdot n', '-\int_A ( J_{F_{GEO}}+J_D + J_{F_{WIND}})', '-J_D','-J_{F_{GEO}}', '-J_{F_{WIND}}', '-J_{D_{SURF}}')
xlabel('Days');
title(['1025.9 < \rho < 1026.3'])
ylabel('m^3/s^{-4}');

% volt = squeeze(sum(sum(sum(mask.*gridvol))));

subplot(2,1,2)
[ax, h1, h2] = plotyy(1:nt, vol, 1:nt, Qm);
grid on
set(get(ax(1), 'YLabel'), 'String', 'Mode Water Volume (m^3)');
set(get(ax(2), 'YLabel'), 'String', 'Q (W m^{-2})');
set(h1, 'LineWidth', 2);
set(h2, 'LineWidth', 2);
set(gcf, 'Color' ,'w')
%%
% 
figure 
plotyy(1:nt, Qta + dJBTAe + dJBTAs + dJFTA + dJFWA + dJENTa, 1:nt, gradient(vol)./vol)
title(num2str(corr(Qta + dJBTAe + dJBTAs + dJFTA + dJFWA + dJENTa, gradient(vol)./vol)))

figure
plot(-dJAzA);
hold on
plot(-dJAxA);
plot(-dJAyA);
plot(Qt)
plot(Qt + dJAxA+dJAyA+dJAzA)
hold off
legend('-dJAzA', '-dJAxA', '-dJAyA', 'Qt', 'SUM');
% hold off
% legend('z', 'x', 'y', 'sum')
% %%
% 
% % plotyy(1:350, Qa , 1:350, JAxA);
% 
% %%
% scatter(QADV, squeeze(nanmean(nanmean(R.*squeeze(mask(:,:,end-1,:))))))

%%
if false
    zeta = OMEGAZ - repmat(f, [1 1 nz nt]);
    cl = [1024 1027];
    for i=1:nt
        subplot(2,1,1)
       pcolor(squeeze(hkpp(:,:,i)).'); shading interp

       colorbar
       hold on
       contour(squeeze(mask(:,:,end-1,i)).', [0 1], 'k');
       hold off
       title(num2str(i));
%        set(gca, 'clim', cl);
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