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
nt = 195;
offset = 0;
% offset = 81; % use for NESEB comparison.
% nt = 35; % use for NESB comparison.
% offset = 0; % Use for NESEB
% nt = 418; % Use for NESEB

forcepath = [pardir 'nesea_frc.nc'];
bdom = false;

% % % GULFZ RUNS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pardir = '/data/thomas/jacob13/GULFZ/';
% basepath = [pardir 'HIS/'];
% path1 = [pardir 'gulfz_grd.nc'];
% forcepath = [pardir 'gulfz_frc.nc'];
% ntp = 5;
% bdom = true;
% % nt =  93;345;
% % offset =  126;

timestring = 'ocean_time';

% % % 500 m high res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = '/groups/thomas2/jacob13/NESEC2/';
ntp = 20;
nt = 30;
offset = 1;
timestring = 'time';
files = dir([basepath,'*his*.nc']);
filespv = dir([basepath,'*pv*.nc']);
    
% Load constants/Initial Params
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
xl =1000:1400;1000:1150;1200:1300;
yl =  300:650;250:350;

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
path = [basepath, files(ceil((0+offset)./ntp)).name];
 
%% PLOTS
if false
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
rho = mask;
dz = mask;
vol = NaN(nt, 1);
Qa = NaN(nt, 1);
negQ = Qa;
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
        JFT= omegazs;
        JBTs= omegazs;
        JBTe = omegazs;
        JFW=omegazs;
        JENT=omegazs;
        tJFN = omegazs;
        tJBN = omegazs;
        Ts = omegazs;
        hkpp = omegazs;
        STRAIN = Ts;
        N2 = omegazs;
        M2 = omegazs;
        Bo = omegazs;
% Qf = NaN(nx, ny, nz, nt);
parfor i=1:nt-offset;
    % 1) Housekeeping
%     disp([num2str(i), '/', num2str((length(files)-ntp)*ntp)]);
    fileind = ceil((i+offset-1)/ntp);
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
    disp(['FileInd: ', num2str(fileind), '   Ind:', num2str(sliceind), '   i=', num2str(i)]);
    

    sliceT = {slice{1}, slice{2}, slice{3},[sliceind sliceind]};

    namepv = filespv(fileind).name;
    pathpv = [basepath, namepv];
    slicePV = {slice{1}, slice{2}, [50 50], sliceT{4}};
   
        [outstruct, outputJV] = ROMSFASTPV_RHS(path,pathpv, pathb, pathf, sliceT, pm, pn, squeeze(sst(:,:,i)),squeeze(sss(:,:,i)), squeeze(Qo(:,:,i)), squeeze(EP(:,:,i)),...
            squeeze(SW(:,:,i)), rho0, Cp, squeeze(tx(:,:,i)), squeeze(ty(:,:,i)), squeeze(tmag(:,:,i)), dxf, dyf,...
            theta_s, theta_b, hc, h, zl, g, f, zmin, ntp, ts, false);
        
        
        mask(:,:,:,i) = outstruct.mask;
%         rho(:,:,:,i) = outstruct.Rho;
        JFT(:,:,i) = outstruct.JFT;
        JBTs(:,:,i) = outstruct.JBTs;
        JBTe(:,:,i) = outstruct.JBTe;
        JFW(:,:,i) = outstruct.JFW;
%         JENT(:,:,i) = outstruct.JENT;
%         Ts(:,:,i ) = squeeze(outstruct.T(:,:,end));
%         STRAIN(:,:,i) = squeeze(outstruct.STRAIN(:,:,end-1));
        hkpp(:,:,i) = outstruct.h;
        vol(i) = outstruct.vol;
        Qat(i) = outstruct.Qat;
        negQ(i) = outstruct.negQ;
        dJAzA(i) = outstruct.dJAzA;
        dJAxA(i) = outstruct.dJAxA;
        dJAyA(i) = outstruct.dJAyA;
        dJBTA(i) = outstruct.dJBTA;
        dJFTA(i) = outstruct.dJFTA;
        dJFWA(i) = outstruct.dJFWA;
        dJBTAs(i) = outstruct.dJBTAs;
        dJBTAe(i) = outstruct.dJBTAe;
        dJBFWe(i) = outstruct.dJBFWe;
%         dJENTa(i) = outstruct.dJENTa;
        hm(i) = outstruct.hm;
        Qm(i) = outstruct.Qm;
%         dz(:,:,:,i) = outstruct.dz;
        Bo(:,:,i) = outstruct.Bo;
        omegazs(:,:,i) = outstruct.omegazs;
%         M2(:,:,i) = outstruct.M2;
%         N2(:,:,i) = outstruct.N2;
%         Qf(:,:,:,i) = outstruct.Q;
%         Tf(:,:,:,i) = outstruct.T;

        %%%% GENERATE J VECTORS
%         outputJV = generateJVectors(pathts, pathuv, sliceTS,outstruct.Bx, outstruct.By,outstruct.Bz, rho0, g,...
%             outstruct.OMEGAX, outstruct.OMEGAY, outstruct.OMEGAZ, outstruct.mask,dxf, dyf,outstruct.dz, -outstruct.alpha, -outstruct.beta);
        
%         outputJV = generateJVectorsRHS(pathpv,path, sliceT, squeeze(outstruct.OMEGAZ),outstruct.OMEGAX, outstruct.OMEGAY,...
%             squeeze(outstruct.Bx(:,:,:)), squeeze(outstruct.By(:,:,:)), 0, 0,...
%             g, squeeze(dxf(:,:,:)), squeeze(dyf(:,:,:)), outstruct.dz, outstruct.mask, ...
%             squeeze(outstruct.Bz(:,:,:)));
        
        dJDn(i) = outputJV.Jda;
        dJFn(i) = outputJV.Jfa;
        
%         dJDt(i) = outputJV.Jda + outputJV.dJBxA + outputJV.dJByA;
%         dJFt(i) = outputJV.Jfa + outputJV.dJFxA + outputJV.dJFyA;
        
        tJFN(:,:,i) = outputJV.JFZ(:,:);
        tJBN(:,:,i) = outputJV.JDZ(:,:);
end
        disp(['Elapsed Time: ', num2str(toc./60)]);

%%
% cl = 1e-13;
% figure
% subplot(5,1,1)
% pcolor(sign(squeeze(nanmean(tJFN, 3))).'.*squeeze(nanmean((tJFN - JFT+2.*JFW).*squeeze(mask(:,:,end-1,:)), 3)).'); shading interp
% set(gca, 'clim', [-1 1].*cl);
% hold on
% contour(squeeze(nanmean(Ts(:,:,:),3).'),10, 'k');
% hold off
% colorbar;
% subplot(5,1,2)
% pcolor(squeeze(nanmean((tJBN - JBTs -1..*JBTe).*squeeze(mask(:,:,end-1,:)), 3)).'); shading interp
% set(gca, 'clim', [-1 1].*cl);
% colorbar;
% subplot(5,1,3)
% pcolor(squeeze(nanmean((JFW).*squeeze(mask(:,:,end-1,:)), 3)).'); shading interp
% set(gca, 'clim', [-1 1].*cl);
% colorbar;
% subplot(5,1,4)
% pcolor(squeeze(nanmean((JENT).*squeeze(mask(:,:,end-1,:)), 3)).'); shading interp
% set(gca, 'clim', [-1 1].*cl);
% colorbar;
% subplot(5,1,5)
% pcolor(squeeze(nanmean((JBTs).*squeeze(mask(:,:,end-1,:)), 3)).'); shading interp
% set(gca, 'clim', [-1 1].*cl);
% colorbar;
%%
% cl = 1e-12;
% ti = 30;
% figure
% subplot(2,1,1)
% pcolor(squeeze((tJFN(:,:,ti) - JFT(:,:,ti))).'); shading interp
% % pcolor(squeeze((JFW(:,:,ti))).'); shading interp
% set(gca, 'clim', [-1 1].*cl);
% hold on
% contour(squeeze(Ts(:,:,ti).'),20, 'k');
% hold off
% colorbar;
% title('$J_F^{NUM} - J_F^{SCALE}$');
% 
% subplot(2,1,2)
% pcolor(squeeze((tJBN(:,:,ti) - 1.2.*JBTs(:,:,ti)-JBTe(:,:,ti))).'); shading interp
% set(gca, 'clim', [-1 1].*cl);
% hold on
% contour(squeeze(Ts(:,:,ti).'), 20, 'k');
% hold off
% colorbar;
% title('$J_D^{NUM} - J_D^{SCALE}$');

%%
figure
subplot(2,1,1)
plot(vol);
title('Mask Volume');
subplot(2,1,2)
yyaxis left
plot(squeeze(nansum(nansum(mask(:,:,end-1,:)))));
yyaxis right
plot(squeeze(nansum(nansum(mask(:,:,2,:)))));
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


QADV = (cumtrapz(dJAzA) + cumtrapz(dJAxA) + cumtrapz(dJAyA)).*ts./vol;

% QADV(m) = cumtrapz(modeltime(m), dJAzA(m) + dJAxA(m) + dJAyA(m));
% QADV = interp1(modeltime(m), QADV(m), modeltime).'./vol;

JBTA = cumtrapz(dJBTA).*ts./vol;
JFTA = cumtrapz(dJFTA).*ts./vol;
JFWA = cumtrapz(dJFWA).*ts./vol;
JBTAs = cumtrapz(dJBTAs).*ts./vol;
JBTAe = cumtrapz(dJBTAe).*ts./vol;
JFWAe = cumtrapz(dJBFWe).*ts./vol;
% JENTA = cumtrapz(dJENTa).*ts./vol;
JDNA = cumtrapz(dJDn).*ts./vol;
JFNA = cumtrapz(dJFn).*ts./vol;

JDZ = (omegazs - repmat(f, [1 1 nt]))./(repmat(f, [1 1 nt])).*JBTs;
dJDZA = squeeze(nansum(nansum(JDZ.*squeeze(mask(:,:,end-1,:)).*repmat(dxf(:,:,end).*dyf(:,:,end), [1 1 nt]))));
JDZA = cumtrapz(dJDZA).*ts./vol;
JDZe = (omegazs - repmat(f, [1 1 nt]))./(repmat(f, [1 1 nt])).*JBTe;
dJDZAe = squeeze(nansum(nansum(JDZe.*squeeze(mask(:,:,end-1,:)).*repmat(dxf(:,:,end).*dyf(:,:,end), [1 1 nt]))));
JDZAe = cumtrapz(dJDZAe).*ts./vol;toc

%% MAKE HISTOGRAM PLOT
% omegavec = reshape(omegazs./repmat(f, [1 1 nt]), [nx.*ny*nt, 1]);
% jbvec = reshape(tJBN.*repmat(dx.*dy, [1 1 nt]), [nx.*ny*nt, 1]);
% jfvec = reshape(tJFN.*repmat(dx.*dy, [1 1 nt]), [nx.*ny*nt, 1]);
% edges = -2:0.5:10;
% [N,edges,bin] = histcounts(omegavec, edges);
% 
% m = bin>0;
% AvgJB = accumarray(bin(m),jbvec(m),[],@nansum);
% AvgJF = accumarray(bin(m),jfvec(m),[],@nansum);
% %% HISTOGRAM FIGURE
% figure
% % subplot(1,2,1);
% bar(0.5.*(edges(2:end)+edges(1:end-1)), -AvgJB)
% ylabel('$-\int J_D^{NUM} \mathrm{d A}$   $\left[m^3s^{-4}\right]$')
% set(gca, 'FontSize', 20, 'XTick', -2:2:10)
% xlabel('$\frac{f+\zeta}{f}$', 'FontSize', 24);
% 
% set(gca, 'xlim', [-2 10]);
% 
% grid on
% set(gcf, 'Color', 'w')

% subplot(1,2,2);
% bar(0.5.*(edges(2:end)+edges(1:end-1)), -AvgJF)
% xlabel('$\frac{f+\zeta}{f}$');
% ylabel('$-\int J_F^{NUM} \mathrm{d A}$')
% set(gca, 'FontSize', 20, 'XTick', -2:2:10)
% set(gca, 'xlim', [-2 10]);
% grid on
% set(gcf, 'Color', 'w', 'Position', [   675   561   848   413])
%%
mt =mod( modeltime./86400, 350);

figure
% subplot(3,1,1)
plot(mt, (Qa+QADV).*vol  , 'LineWidth', 2);
hold on
% plot(smooth((Qa +QADV).*vol, 2), 'LineWidth', 2);

plot(mt, -JFTA.*vol, 'LineWidth', 2);
plot(mt, -(JBTAe+JBTAs+0*JDZA).*vol, 'LineWidth', 2);
plot(mt, -JFWA.*vol, 'LineWidth', 2)


plot(mt, -(JFTA+JBTAe+1.2*JBTAs+JFWA+0*JDZAe+0*JDZA).*vol, 'LineWidth', 2, 'LineStyle', '--')
% plot(modeltime, -JENTA.*vol, '--');
plot(mt, -(JDNA).*vol, '--g');
plot(mt, -JFNA.*vol, '--b');
plot(mt, -(JDNA+JFNA).*vol, '--r', 'LineWidth', 2);
% plot(-(JFzA + JBzA+JFxA + JFyA), 'k','LineWidth', 2)
hold off
legend('$\Delta \int_V Q + \int_t\int_\Sigma uQ \cdot n$','$-\int_t\int_A J_{F_{GEO}}$',...
    '$-\int_t\int_A J_{D}$', '$-\int_t\int_A J_{F_{WIND}}$', '$-\int_t\int_A (J_{F_{GEO}} + J_D + J_{F_{WIND}})$',...
    '$-\int_t\int_A J_D^{NUM}$', '$-\int_t\int_A J_F^{NUM}$', '$-\int_t\int_A(J_D^{NUM}+J_F^{NUM})$', 'Location', 'NorthWest')
grid on
title(['$1025.9 < \rho < 1026.3$'])

ylabel('$m^3 s^{-3}$');
set(gca, 'FontSize', 16);
% datetick('x')
% set(gca, 'xlim', [mt(1) mt(35)]);
xlabel('Yearday');
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
set(gcf, 'Color', 'w')

%% NON-CONS d/dt terms.
s =1;
Qta = gradient((Qa+QADV).*vol, ts);

Qtana = gradient(Qa.*vol, ts);
% Qta = smooth(Qt + dQADV, s);
figure
% subplot(2,1,1)
plot(mt, smooth(Qta,s), 'LineWidth', 2);
hold on
plot(mt, -smooth(dJBTAe  + 1.2*dJBTAs +dJFTA+1.*dJFWA+0*dJDZAe+0*dJDZA,s), 'LineWidth', 2);
plot(mt, -smooth(dJBTA, s));
plot(mt, -smooth(dJFTA,s));
plot(mt, -smooth(dJFWA,s), 'r');
plot(mt, -smooth(dJBTAs,s),'--')
plot(mt, -smooth(dJENTa,s))
% plot(mt, -smooth(dJDZAe,s))

plot(mt, -smooth(dJDn, s));
plot(mt, -smooth(dJFn, s));
plot(mt, -smooth(dJFn+dJDn, s), 'k')
plot(mt, smooth(Qtana, s), 'cy')
% df = 0.*dJBTAe + 1.*dJENTa + dJBTAs + 2/3.*(-dJFWA-dJFTA);
% ff = dJFTA + dJFWA;
% plot(-smooth(df,s));
% plot(-smooth(ff,s));
% plot(df + ff)
% plot(-smooth(dJBTA+dJFTA,s));
hold off
grid on
legend('$\partial/\partial t \int_V Q + \int_\Sigma uQ \cdot n$', '$-\int_A ( J_{F_{GEO}}+J_D + J_{F_{WIND}})$',...
    '$-J_D$','$-J_{F_{GEO}}$', '$-J_{F_{WIND}}$', '$-J_{D_{SURF}}$', '$J_{D_{ENT}}$', '$J_D^N$', '$J_F^N$', 'NUMSUM', 'location', 'NorthEastOutside')
xlabel('Days');
title(['$1025.9 < \rho < 1026.3$'])
ylabel('$m^3/s^{-4}$');
set(gca, 'xlim', [mt(1) mt(end)]);

set(gcf, 'Color', 'w')
% volt = squeeze(sum(sum(sum(mask.*gridvol))));

% subplot(2,1,2)
% [ax, h1, h2] = plotyy(1:nt, vol, 1:nt, Qm);
% grid on
% set(get(ax(1), 'YLabel'), 'String', 'Mode Water Volume (m^3)');
% set(get(ax(2), 'YLabel'), 'String', 'Q (W m^{-2})');
% set(h1, 'LineWidth', 2);
% set(h2, 'LineWidth', 2);
% set(gcf, 'Color' ,'w')
%%
figure
subplot(1,2,1)
scatter(dJFn, dJFTA);
onetoone
grid on
cr = corr(dJFn, dJFTA, 'rows', 'pairwise');
title(num2str(cr))
subplot(1,2,2)
scatter(dJDn, dJBTAe + 1.2*dJBTAs)
cr = corr(dJDn(5:end), dJBTAe(5:end) + 1.2*dJBTAs(5:end)+0*(dJDZA(5:end)+dJDZAe(5:end)), 'rows', 'pairwise');
title(num2str(cr))
onetoone
grid on

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


if false
%%
lat = ncread('/groups/thomas1/jacob13/GulfStream/NESEA/nesea_grd.nc', 'lat_rho');
lon = ncread('/groups/thomas1/jacob13/GulfStream/NESEA/nesea_grd.nc', 'lon_rho');

%% RO T plots
ti=4;
tlim = [3 25];

Ro = (squeeze(omegazs(:,:,ti)) - f)./f;
Rib = abs(f.^2.*squeeze(N2(:,:,ti))./squeeze(M2(:,:,ti)));
h  = ones(3); h(2,2) = 2; h = h./sum(sum(h));
% Rib = filter2(h, Rib);
Ro = imgaussfilt(Ro, [1  1]);
Rib = 1./imgaussfilt(Rib, [2 2]);
subplot(1,2,1)
pcolor( (0:2001)./2,(0:1601)./2, squeeze(Ts(:,:,ti)).'); shading interp
hold on
% contour((0:2001)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
hold off
set(gca, 'clim', tlim);
% title('$Temperature$');
set(gca, 'FontSize', 16);
axis equal
xlabel('km'); ylabel('km');
set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 
% colormap(gca, cptcmap('BlueWhiteOrangeRed.cpt'))
set(gcf, 'Color', 'w', 'Position', [   675    63   841   911]);

subplot(1,2,2)
pcolor( (0:2001)./2,(0:1601)./2, squeeze(Ro(:,:)).'); shading interp
hold on
% contour((0:2001)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
hold off
set(gca, 'clim', [-1 1]*2);
% title('$Temperature$');
set(gca, 'FontSize', 16);
axis equal
xlabel('km'); ylabel('km');
set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 
colormap( cptcmap('BlueWhiteOrangeRed.cpt'))
set(gcf, 'Color', 'w', 'Position', [   675    63   841   911]);


%% Ratio EBF Plots
ti=4;
tlim = [3 25];
clr = [0 5];
EBFgtoW = 0.1.*squeeze(hkpp(:,:,ti).*sqrt(M2(:,:,ti)))./(f.*sqrt(squeeze(tmag(:,:,ti))./1027));
EBFgtoB = 0.05.*squeeze(hkpp(:,:,ti).^2.*M2(:,:,ti))./(f.*squeeze(abs(Bo(:,:,ti))));
h  = ones(3); h(2,2) = 2; h = h./sum(sum(h));
% Rib = filter2(h, Rib);
% Ro = imgaussfilt(Ro, [1  1]);
% Rib = 1./imgaussfilt(Rib, [2 2]);
figure
% subplot(1,2,1)
pcolor( (0:2001)./2,(0:1601)./2, squeeze(EBFgtoW).'); shading interp
hold on
% contour((0:2001)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
hold off
set(gca, 'clim', clr);
% title('$Temperature$');
set(gca, 'FontSize', 18);
axis equal
xlabel('km'); ylabel('km');
set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 
% colormap(gca, cptcmap('BlueWhiteOrangeRed.cpt'))
cb = colorbar('SouthOutside');
colormap( cptcmap('haline.cpt'))
set(gcf, 'Color', 'w', 'Position', [  170          21        1309         944]);
t= text(50, 750,0, '$\frac{\tau_g}{\tau_w}$', 'FontSize', 40, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');
set(cb, 'Ticks', 0:1:5);
set(cb, 'FontSize', 24)
%%
title('$\frac{EBF_g}{EBF} \sim 0.1 \frac{H |\nabla_h b|}{f \sqrt{\tau/\rho}}$', 'Interpreter', 'Latex', 'FontSize', 20);
subplot(1,2,2)
pcolor( (0:2001)./2,(0:1601)./2, squeeze(EBFgtoB(:,:)).'); shading interp
hold on
% contour((0:2001)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
hold off
set(gca, 'clim', clr);
% title('$Temperature$');
set(gca, 'FontSize', 16);
axis equal
xlabel('km'); ylabel('km');
set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 
colormap( cptcmap('haline.cpt'))

title('$\frac{EBF_g}{Bo} \sim 0.05 \frac{H^2 |\nabla_h b|^2}{f|B_o|}$', 'Interpreter', 'Latex', 'FontSize', 20);
cb = colorbar;
set(gcf, 'Color', 'w', 'Position', [  170          21        1309         944]);
set(cb, 'Position',[0.9253    0.3305    0.0213    0.3741]);
% subplot(1,2,2)
% pcolor( (0:2001)./2,(0:1601)./2, squeeze(Rib(:,:)).'); shading interp
% hold on
% % contour((0:2001)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
% hold off
% set(gca, 'clim', [-1 1]*2);
% % title('$Temperature$');
% set(gca, 'FontSize', 16);
% axis equal
% xlabel('km'); ylabel('km');
% set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 
% colormap(gca, cptcmap('BlueWhiteOrangeRed.cpt'))
% set(gcf, 'Color', 'w', 'Position', [   675    63   841   911]);
%% B Gradient ONLY
Bgradn = squeeze(sqrt(M2(:,:,ti)))./f.^2;
pcolor( (0:2001)./2,(0:1601)./2, squeeze(Bgradn).'); shading interp
hold on
% contour((0:2001)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
hold off
set(gca, 'clim',[0 100]);
% title('$Temperature$');
set(gca, 'FontSize', 16);
axis equal
xlabel('km'); ylabel('km');
set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 
colormap( cptcmap('haline.cpt'))

title('$\frac{|\nabla_h b|}{f^2}$', 'Interpreter', 'Latex', 'FontSize', 26);
cb = colorbar;
set(gcf, 'Color', 'w', 'Position', [           170         245        1112         719]);
set(cb, 'Position',[       0.8860    0.1090    0.0286    0.8076]);
set(cb, 'Ticks', [0 100]);
%% TEMPERATURE ONLY
pcolor( (0:2001)./2,(0:1601)./2, squeeze(Ts(:,:,ti)).'); shading interp
hold on
% contour((0:2001)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
hold off
set(gca, 'clim', tlim);
% title('$Temperature$');
set(gca, 'FontSize', 16);
axis equal
xlabel('km'); ylabel('km');
set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 
colormap(gca, cptcmap('BlueWhiteOrangeRed.cpt'))
set(gcf, 'Color', 'w', 'Position', [   675    63   841   911]);

%%
jstring = 'BlWhRe.cpt';
ti = 3;
cl = 2e-12;
mc = .9; %coast color
tlim = [3 25];
gap = [.12 .2]; margh = .1; margw=.1;
figure
% subtightplot(2, 2,1, gap, margh, margw);
% pcolor( (0:2001)./2,(0:1601)./2, squeeze(Ts(:,:,ti)).'); shading interp
% hold on
% % contour((0:2001)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
% hold off
% set(gca, 'clim', tlim);
% title('$Temperature$');
% set(gca, 'FontSize', 16);
% axis equal
% xlabel('km'); ylabel('km');
% set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 
% colormap(gca, cptcmap('BlueWhiteOrangeRed.cpt'))

subtightplot(1, 2,1, gap, margh, margw);
pcolor( (0:2001)./2,(0:1601)./2, -squeeze(tJFN(:,:,ti)).'); shading interp
hold on
% contour((0:2001)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
hold off
set(gca, 'clim', [-1 1].*cl);
title('$-J_F^{NUM}$');
set(gca, 'FontSize', 16);
axis equal
xlabel('km'); ylabel('km');
set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 
colormap(gca, cptcmap(jstring));

subtightplot(1, 2,2, gap, margh, margw);
pcolor( (1:2002)./2,(1:1602)./2, -squeeze(tJBN(:,:,ti)).'); shading interp
hold on
% contour((1:2002)./2, (1:1602)./2, squeeze(Ts(:,:,ti)).', 0:2:30, 'k')
hold off
set(gca, 'clim', [-1 1].*cl);
title('$-J_D^{NUM}$');
set(gca, 'FontSize', 16);
axis equal
xlabel('km'); ylabel('km');
cb = colorbar;
set(get(cb, 'YLabel'), 'String', '$ms^{-4}$', 'Interpreter', 'Latex');
set(gca, 'xlim', [0 2001]./2, 'ylim', [0 1601]./2); 

% subtightplot(2, 1,2, gap, margh, margw);
% axesm('mercator','MapLatLimit', [34 44], 'MapLonLimit', [-72.5 -58], ...
%     'ParallelLabel', 'on', 'MeridianLabel', 'on','PlabelLocation', 5, 'MLabelLocation', 5,'MLabelParallel', 'South',...
%     'PLineLocation', 5, 'MLineLocation', 5, 'Frame', 'on');
% % axesm('mercator')
% framem; gridm; axis on; tightmap;
% pcolorm( squeeze(lat(:,:)),squeeze(lon(:,:)), -squeeze(tJBN(:,:,ti))); shading interp
% contourm(lat, lon, squeeze(Ts(:,:,ti)), 0:2:30, 'k')
% hold off
% % geoshow(shorelines(land), 'DisplayType' , 'Polygon', 'FaceColor', [1 1 1].*mc)
% set(gca, 'clim', [-1 1].*cl);
% t= textm(45, 360-82.5,0, '$T$', 'FontSize', 20, 'Color', 'k', ...
%     'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');
% title('$-J_D^{NUM}$')
% colormap(cptcmap('dkbluered.cpt'));
colormap(gca, cptcmap(jstring));

set(gcf, 'Color', 'w', 'Position', [   675    63   841   911]);
set(cb, 'Position', [  0.7589    0.3233    0.0314    0.3403]);

end