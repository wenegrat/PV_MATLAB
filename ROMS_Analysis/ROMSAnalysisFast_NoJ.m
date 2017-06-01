% NOTING all things to be looked at later:
% XXX - grid is spherical.
% XXX - Probably should directly calculate Pressure gradients (not press)
tic

%%
% CAN REPRODUCE THE ISSUE IN NESEA using nt = 43; offset = 150; xl =
% 800:1050; yl = 300:500. Working theory is advection of GS meander into
% domain is not captured by fluxes at edge of isopycnal.

%PARAMETERS TO CHANGE'

% % % 500 m run - initial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pardir = '/groups/thomas1/jacob13/GulfStream/NESEA/';
basepath = [pardir 'HIS/'];
path1 = '/groups/thomas1/jacob13/GulfStream/NESEA/HIS/nesea_his.0000.nc';
ntp = 2;
offset = 10;
nt = 180;190;

% offset = 81; % use for NESEB comparison.
% nt = 35; % use for NESB comparison.
% offset = 0; % Use for NESEB
% nt = 418; % Use for NESEB

forcepath = [pardir 'nesea_frc.nc'];
bdom = false;

% % % GULFZ RUNS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pardir = '/data/thomas/jacob13/GULFZ/';
basepath = [pardir 'HIS/'];
path1 = [pardir 'gulfz_grd.nc'];
forcepath = [pardir 'gulfz_frc.nc'];
ntp = 5;
bdom = false;
nt =  240;
offset =  60;
% 
timestring = 'ocean_time';
faststring = false;

% % % 500 m high res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XXXXX----- USE ROMSAnalysisFast.m
% basepath = '/groups/thomas1/jacob13/NESEC2/';
% bdom = false;
% ntp = 10;
% nt = 50;
% offset = 10;
% timestring = 'time';
% filesuv = dir([basepath,'*uv*.nc']);
% filests = dir([basepath,'*ts*.nc']);
% faststring = true;

% Load constants/Initial Params
files = dir([basepath,'*his*.nc']);
path = [basepath, files(1).name];
time = ncread(path, timestring); %XXXX
% nt = length(time);
ts = time(2)-time(1);


% startmodel = startdate + datenum(0, 0, 0, 0, 0, tf(1)-oti);
% dateoff = datenum(0, 0, offset);
% modeltime = inidate + dateoff + datenum(0, 0, 0:1:nt-1);


% modeltime = datenum(2012,8,(26+offset):1:(26+offset+(nt)-1), 0,0,0);


%

zl = 1:50;
nz = length(zl);
xl =  1:1602;1:2002;100:300;1000:1150;1200:1300;
yl =   1:922;1:802; 1:150;250:350;

% USE THESE FOR LARGE DOMAIN PV BUDGET
xl = 600:1602;
yl = 1:550;

% USE THIS FOR LARGE DOMAIN TO OVERLAP.
if bdom
[xl, yl, offset, nt] = findlimitsLarge(xl, yl, nt);
end

if offset>0
tf = ncread([basepath, files(ceil((offset+1)./ntp)).name], timestring); %XXXX
else
    tf = ncread([basepath, files(1).name], timestring); % XXXXXX
end
modeltime = tf(1):ts:((nt-1)*ts+tf(1));

% xl = 1237:1336;
% yl = 295:395;
slice =  {[xl(1) xl(end)], [yl(1) yl(end)], [zl(1) zl(end)], 0};



path = path1;
f = ncread(path, 'f');
[nx ny] = size(f)
f = f(xl, yl);
rho0 = 1027.4;
g = 9.81;
Cp = 3985; % Value per Jon's email 4/20/17
zmin = 5.86;
pm = ncread(path, 'pm'); %XXXXXXX
pm = pm(xl, yl);
% pmmetric = repmat(pm, [1 1 nz nt]);

pn = ncread(path, 'pn'); % XXXXXXX
pn = pn(xl,yl);
dx = 1./pm;
dy = 1./pn;
dxf = repmat(dx, [1 1 nz]);
dyf = repmat(dy, [1 1 nz]);

h = ncread(path, 'h'); %XXXXXXXXX

latst = ncread(path1, 'lat_rho'); %XXXXXXX
lats = latst(xl, yl);
lonst = ncread(path1, 'lon_rho'); %XXXXXXXXXXX
lons = lonst(xl, yl);
disp(['Sub-domain borders (clockwise from bottom left): (', num2str(lats(1,1)), ',',num2str(lons(1,1)),')  (',...
    num2str(lats(1,end)),',',num2str(lons(1,end)),')  (', num2str(lats(end,end)),',',num2str(lons(end,end)),')  (',...
    num2str(lats(end,1)),',',num2str(lons(end,1)),')']);
% FULL DOMAIN
path = [basepath, files(ceil((1+offset)./ntp)).name];
 
%% PLOTS
if true
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
% axesm('mercator', 'MapLatLimit', [26.75 47], 'MapLonLimit', [-90 -55],...
%     'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5, 'PLineLocation', 5, 'MLineLocation', 5);
% % axesm('mercator')
% framem; gridm; axis on; tightmap;
% pcolorm( squeeze(latst(:,:)),squeeze(lonst(:,:)),(squeeze(T(:,:,end,1)))); shading interp
% 
lonf = [lons(1,:) lons(:,end).' fliplr(lons(end,:)) flipud(lons(:,1)).'];
latf = [lats(1,:) lats(:,end).' fliplr(lats(end,:)) flipud(lats(:,1)).'];
% hold on
% hf = fillm(latf, lonf, 'k');
% set(hf, 'FaceColor', 'none');
% set(hf, 'LineWidth', 2);
% % set(h, 
% % rectangle('Position', [xl(1) yl(1) xl(end)-xl(1) yl(end)-yl(1)]);
% hold off

%LOAD SHORELINE
shorelines = gshhs('~/GSHH/gshhs_c.b');%_h for high, _c for crude, _i for intermediate, from https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
levels = [shorelines.Level];
land = (levels==1);

cl = 5e-12;
mc = .9; %coast color
gap = [.05 .01]; margh = .1; margw=.1;
figure
axesm('mercator','MapLatLimit', [26.75 47], 'MapLonLimit', [-85 -56.5], ...
    'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5,'MLabelParallel', 'South',...
    'PLineLocation', 5, 'MLineLocation', 5, 'Frame', 'on');
% axesm('mercator')
framem; gridm; axis on; tightmap;
pcolorm( squeeze(latst(:,:)),squeeze(lonst(:,:)),(squeeze(T(:,:,end,1)))); shading interp
hold on
% contourm(latst, lonst, squeeze(T(:,:,end,1)), 0:1:30, 'k')
hf = fillm(latf, lonf, 'k');
set(hf, 'FaceColor', 'none');
set(hf, 'LineWidth', 2);
hold off
geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc)
set(gca, 'clim', [10 28]);
t= textm(45, 360-82.5,0, '$T$', 'FontSize', 20, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');

% colorbar;
end

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
SSS = ncread(forcepath, 'SSS');
SSS = SSS(xl, yl,:);
Qor = repmat(Qor, [1 1 2]);
EPr = repmat(EPr, [1 1 2]);
tau_u = repmat(tau_u, [1 1 2]);
tau_v = repmat(tau_v, [1 1 2]);
SST = repmat(SST, [1 1 2]);
SSS = repmat(SSS, [1 1 2]);
SWR = repmat(SWR, [1 1 2]);
% climtime = datenum(2012,1, 15:30:(360*2));
ds = 86400;
climtime = ((360+15)*ds):30*ds:(360*2*ds+360*ds);
% disp('here')
[nx, ny, ~] = size(Qor);

% Qo = NaN(nx, ny, nt);
% EP = Qo;
% SW = Qo;
% tx = Qo;
% sst = Qo;
% sss = Qo;
% ty = Qo;

[X, Y, T] = ndgrid(xl, yl, climtime);
[Xm, Ym, mt] = ndgrid(xl, yl, modeltime);
QG = griddedInterpolant(X, Y, T, Qor, 'cubic', 'none');
EPG = griddedInterpolant(X, Y, T, EPr, 'cubic', 'none');
TXG = griddedInterpolant(X, Y, T, tau_u, 'cubic', 'none');
TYG = griddedInterpolant(X, Y, T, tau_v, 'cubic', 'none');
SSTG = griddedInterpolant(X, Y, T, SST,'cubic', 'none');
SWG = griddedInterpolant(X, Y, T, SWR, 'cubic', 'none');
SSSG = griddedInterpolant(X, Y, T, SSS, 'cubic', 'none');

Qo = QG(Xm, Ym, mt);
EP = EPG(Xm, Ym, mt);
tx = TXG(Xm, Ym, mt);
ty = TYG(Xm, Ym, mt);
sst = SSTG(Xm, Ym, mt);
SW = SWG(Xm, Ym, mt);
sss = SSSG(Xm, Ym, mt);

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
dJFn = Qa;
dJDn = Qa;
dJFt = Qa;
dJDt = Qa;

hm = Qa;
Qm = hm;
% Bo = hm;
omegazs = NaN(nx, ny, nt);
omegaf = NaN(nx, ny, nz, nt);
        JFT= omegazs;
        JBTs= omegazs;
        JBTe = omegazs;
        JFW=omegazs;
        JENT=omegazs;
%         tJFN = omegazs;
%         tJBN = omegazs;
        Ts = omegazs;
        hkpp = omegazs;
        bo = omegazs;
        STRAIN = Ts;
        M2 = Ts;
% Qf = NaN(nx, ny, nz, nt);

parfor i=1:nt;
    % 1) Housekeeping
    disp([num2str(i), '/', num2str(nt)]);

     if offset>0
          fileind = ceil((i+offset)/ntp);
     else
          fileind = ceil(i./ntp)
     end

    name = files(fileind).name;
    path = [basepath, name];
     
    if fileind==1;
        fileindb = 1;
    else
        fileindb = fileind-1;
    end
    pathb = [basepath, files(fileindb).name];
    
    if fileind==length(files)
      fileindf = fileind;
    else
        fileindf =  fileind+1;
    end
    pathf = [basepath, files(fileindf).name];
    
    sliceind = mod(i,ntp);
    if sliceind==0; 
        sliceind=ntp;
    end
%     fileindts = fileind-1;
%     sliceindts = sliceind-1;
%     if sliceind==1;
%        fileindts =  fileind-2;
%        sliceindts = ntp;
%     end
%     if (i==1)
%         fileindts = fileind-1;
%         sliceindts = 1;
%     end
    sliceT = {slice{1}, slice{2}, slice{3},[sliceind sliceind]};

%     nameuv = filesuv(fileindts).name;
%     pathuv = [basepath, nameuv];
%     namets = filests(fileindts).name;
%     pathts = [basepath, namets];
%     sliceTS = {slice{1}, slice{2}, slice{3},[sliceindts sliceindts]};
   
        outstruct = ROMSFASTPV(path, pathb,pathf, sliceT, pm, pn, squeeze(sst(:,:,i)),squeeze(sss(:,:,i)), squeeze(Qo(:,:,i)), squeeze(EP(:,:,i)),...
            squeeze(SW(:,:,i)), rho0, Cp, squeeze(tx(:,:,i)), squeeze(ty(:,:,i)), squeeze(tmag(:,:,i)), dxf, dyf,...
            theta_s, theta_b, hc, h, zl, g, f, zmin, ntp, ts, true);
        
        
        mask(:,:,:,i) = outstruct.mask;
        JFT(:,:,i) = outstruct.JFT;
        JBTs(:,:,i) = outstruct.JBTs;
        JBTe(:,:,i) = outstruct.JBTe;
        JFW(:,:,i) = outstruct.JFW;
        JENT(:,:,i) = outstruct.JENT;
        Ts(:,:,i ) = squeeze(outstruct.T(:,:,end));
        STRAIN(:,:,i) = squeeze(outstruct.STRAIN(:,:,end-1));
        M2(:,:,i) = squeeze(outstruct.M2(:,:));
        hkpp(:,:,i) = outstruct.h;
        bo(:,:,i) = outstruct.Bo;
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
        omegaf(:,:,:,i) = outstruct.OMEGAZ;
%         Qf(:,:,:,i) = outstruct.Q;
%         Tf(:,:,:,i) = outstruct.T;

        %%%% GENERATE J VECTORS
%         outputJV = generateJVectors(pathts, pathuv, sliceTS,outstruct.Bx, outstruct.By,outstruct.Bz, rho0, g,...
%             outstruct.OMEGAX, outstruct.OMEGAY, outstruct.OMEGAZ, outstruct.mask,dxf, dyf,outstruct.dz, -outstruct.alpha, -outstruct.beta);
%         dJDn(i) = outputJV.Jda;
%         dJFn(i) = outputJV.Jfa;
%         
%         dJDt(i) = outputJV.Jda + outputJV.dJBxA + outputJV.dJByA;
%         dJFt(i) = outputJV.Jfa + outputJV.dJFxA + outputJV.dJFyA;
%         
%         tJFN(:,:,i) = outputJV.Jf(:,:,end-1);
%         tJBN(:,:,i) = outputJV.Jd(:,:,end-1);
end
        disp(['Elapsed Time: ', num2str(toc./60)]);

%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do post calculations
% vol = nanmean(vol);
Qt = gradient(Qat, ts); %Take the time derivative

% Try to work around volume issues.
% m = smooth(abs(gradient(vol)./vol),1) < 0.05;
% Qa = NaN(length(Qt),1);
% Qa(m) = cumtrapz(modeltime(m).',Qt(m));
% Qa = interp1(modeltime(m), Qa(m), modeltime).'./vol;

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

JDZ = (omegazs )./(repmat(f, [1 1 nt])).*JBTs;
dJDZA = squeeze(nansum(nansum(JDZ.*squeeze(mask(:,:,end-1,:)).*repmat(dxf(:,:,end).*dyf(:,:,end), [1 1 nt]))));
JDZA = cumtrapz(dJDZA).*ts./vol;

JBZ = (omegazs )./(repmat(f, [1 1 nt])).*JBTe;
dJBZA = squeeze(nansum(nansum(JBZ.*squeeze(mask(:,:,end-1,:)).*repmat(dxf(:,:,end).*dyf(:,:,end), [1 1 nt]))));
JBZA = cumtrapz(dJBZA).*ts./vol;

JFZ = (omegazs - repmat(f, [1 1 nt]))./(repmat(f, [1 1 nt])).*JFT;
dJFZA = squeeze(nansum(nansum(JFZ.*squeeze(mask(:,:,end-1,:)).*repmat(dxf(:,:,end).*dyf(:,:,end), [1 1 nt]))));
JFZA = cumtrapz(dJFZA).*ts./vol;
toc
%%

%%
figure
% subplot(3,1,1)
plot(modeltime, (Qa+QADV).*vol  , 'LineWidth', 2);
hold on
% plot(smooth((Qa +QADV).*vol, 2), 'LineWidth', 2);

plot(modeltime, -JFTA.*vol, 'LineWidth', 2);
plot(modeltime, -(JBTAe+JBTAs+0.*JENTA).*vol, 'LineWidth', 2);
plot(modeltime, -JFWA.*vol, 'LineWidth', 2)



plot(modeltime, -(JFTA+(JBTAe+1.2.*JBTAs)+0.*(1.2.*JDZA+JBZA)*.1./.15+JFWA).*vol, 'LineWidth', 2, 'LineStyle', '--')
% plot(modeltime, -JENTA.*vol, '--');
% plot(modeltime, -(JDNA).*vol, '--g');
% plot(modeltime, -JFNA.*vol, '--b');
% plot(modeltime, -(JDNA+JFNA).*vol, 'r');
% plot(-(JFzA + JBzA+JFxA + JFyA), 'k','LineWidth', 2)
hold off
legend('$\Delta \int_V Q + \int_t\int_\Sigma uQ \cdot n$','$-\int_t\int_A J_{F_{GEO}}$', '$-\int_t\int_A J_{D}$', '$-\int_t\int_A J_{F_{WIND}}$', '$-\int_t\int_A (J_{F_{GEO}} + J_D + J_{F_{WIND}})$','Location', 'NorthWest')
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

%% NON-CONS d/dt terms ALL
s =1;
Qta = gradient((Qa+1.*QADV).*vol, ts);

Qtana = gradient(Qa.*vol, ts);
% Qta = smooth(Qt + dQADV, s);
figure
% subplot(2,1,1)
plot(smooth(Qta,s), 'LineWidth', 2);
hold on
plot(-smooth(1.*(dJBTAe  + 1.2*dJBTAs) +dJFTA+1.*dJFWA+0.*dJENTa+0.*(1.*dJBZA+1.2*dJDZA).*0.1./0.15,s), 'LineWidth', 2);
plot(-smooth(dJBTA, s));
plot(-smooth(dJFTA,s));
plot(-smooth(dJFWA,s), 'r');
plot(-smooth(dJBTAs,s),'--')
plot(-smooth(dJBZA,s))
plot(-smooth(dJDZA, s), 'g');

hold off
grid on
legend('$\partial/\partial t \int_V Q + \int_\Sigma uQ \cdot n$', '$-\int_A ( J_{F_{GEO}}+J_D + J_{F_{WIND}})$',...
    '$-J_D$','$-J_{F_{GEO}}$', '$-J_{F_{WIND}}$', '$-J_{D_{SURF}}$', '$J_{D_{ENT}}$', '$J_D^N$', '$J_F^N$', 'NUMSUM')
xlabel('Days');
title(['1025.9 < \rho < 1026.3'])
ylabel('m^3/s^{-4}');

%% NON-CONS d/dt terms FANCY
gap = [.1 .1]; margh = .1; margw = .12;
s =5;
Qta = gradient((Qa+1.*QADV).*vol, ts);
Qtana = gradient(Qa.*vol, ts);
% Qta = smooth(Qt + dQADV, s);
mn = modeltime./(86400) + mod(modeltime(1)./(86400), 360) - modeltime(1)./86400;
mn = mn./30;
% mn = 
figure
subtightplot(2,1,1, gap, margh, margw)
plot(mn, smooth(Qta,s), 'LineWidth', 3);
hold on
plot(mn, -smooth(dJBTA, s), 'LineWidth', 2);
plot(mn, -smooth(dJFTA,s), 'LineWidth', 2);
plot(mn, -smooth(dJFWA,s), 'LineWidth', 2);
plot(mn, -smooth(dJBTAs,s),'--', 'LineWidth', 2);
plot(mn, -smooth(1.*(dJBTAe  + 1.2*dJBTAs) +dJFTA+1.*dJFWA,s), 'LineWidth', 3);

% plot(-smooth(dJBZA,s))
% plot(-smooth(dJDZA, s), 'g');

hold off
grid on
 legend('$\partial/\partial t \int_V Q + \int_\Sigma uQ \cdot n$', '$-\int_A ( J_{F_{GEO}}+J_D + J_{F_{WIND}})$',...
     '$-J_D$','$-J_{F_{GEO}}$', '$-J_{F_{WIND}}$', '$-J_{D_{SURF}}$', '$J_{D_{ENT}}$', '$J_D^N$', '$J_F^N$', 'NUMSUM')
xlabel('Month');
% title(['$1025.9 < \rho < 1026.3$'])
title('Mode Water PV Budget')
ylabel('$m^3s^{-4}$');
set(gca, 'FontSize', 18)

t = text(10.25, .15, 'Rate of Change');
set(t, 'EdgeColor', 'k', 'FontSize', 16, 'BackgroundColor','w')
set(gca, 'xlim', [10 mn(end)], 'XTick', 10:1:18, 'XTickLabel', {10, 11, 12, 1, 2, 3, 4, 5, 6, 7});

subtightplot(2,1,2, gap, margh, margw)
plot(mn, (Qa+QADV).*vol  , 'LineWidth', 3);
hold on
% plot(mn, -(JBTA).*vol, 'LineWidth', 2);
% plot(mn, -JFTA.*vol, 'LineWidth', 2);
% plot(mn, -JFWA.*vol, 'LineWidth', 2)
% plot(mn, -(JBTAe + JFTA + JFWA+1.2.*JBTAs).*vol, 'LineWidth', 2);
% plot(mn, (Qa+QADV+1.2.*JBTAs).*vol  , 'LineWidth', 2);
set(gca, 'ColorOrderIndex', 5)
plot(mn, -1.2.*JBTAs.*vol, '--','LineWidth', 2);
hold off
grid on
% title(['1025.9 < \rho < 1026.3'])
% title('Cumulative PV Budget')
t = text(10.25, 1e5, 'Cumulative $\Delta q$');
set(t, 'EdgeColor', 'k', 'FontSize', 16, 'BackgroundColor','w')
set(gca, 'xlim', [10 mn(end)], 'XTick', 10:1:18, 'XTickLabel', {10, 11, 12, 1, 2, 3, 4, 5, 6, 7});
xlabel('Month');

ylabel('$m^3 s^{-3}$');
set(gca, 'FontSize', 18);
set(gcf, 'Color', 'w', 'Position', [675   318   757   656]);

%%
% % 
% figure 
% plotyy(1:nt, Qta + dJBTAe + dJBTAs + dJFTA + dJFWA + dJENTa, 1:nt, gradient(vol)./vol)
% title(num2str(corr(Qta + dJBTAe + dJBTAs + dJFTA + dJFWA + dJENTa, gradient(vol)./vol)))
% 
% figure
% plot(-dJAzA);
% hold on
% plot(-dJAxA);
% plot(-dJAyA);
% plot(Qt)
% plot(Qt + dJAxA+dJAyA+dJAzA)
% hold off
% legend('-dJAzA', '-dJAxA', '-dJAyA', 'Qt', 'SUM');
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
    zeta = OMEGAZ - repmat(f, [1 1 nt]);
    cl = [1024 1027];
    for i=1:nt
        subplot(2,1,1)
%        pcolor(squeeze(hkpp(:,:,i)).'); shading interp

       colorbar
       hold on
       contour(squeeze(mask(:,:,end-1,i)).', [0 1], 'k');
       hold off
       title(num2str(i));
%        set(gca, 'clim', cl);
       subplot(2,1,2)
       pcolor((squeeze(zeta(:,:,i))./f).'); shading interp
          hold on
       contour(squeeze(mask(:,:,i)).', [0 1], 'k');
       hold off
       set(gca, 'clim', [-2 2]);
       colorbar
       pause(0.001)
    end
end