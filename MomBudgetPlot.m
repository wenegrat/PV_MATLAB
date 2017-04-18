%% MOMENTUM BUDGET PLOT
ti = 188; 
te = 196;
% ti = 189;
% te = 193;
yi = 181;126; %62.75 km
xi = 105;  %41.75 km

xi = 21;
yi=120;
% xi = 86;
nx = 21; ny =21; nt = te-ti+1;
ind = ceil(nx/2);
slice = {[xi xi], [yi yi], 0, [ti te]};
sliceD = {[xi-floor(nx/2) xi+floor(nx/2)], [yi-floor(ny/2) yi+floor(ny/2)], 0, [ti te]};

sliceEta = {sliceD{1}, sliceD{2}, [1 1], sliceD{4}};


ztmp = ncread(statefile, 'Z');
ztmp = ztmp(1:nz);
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);
zL = length(ztmp);

% Tendency Term
dudt = squeeze(GetVar(statefile, diagfile, {'TOTUTEND', '(1)./86400'}, sliceD));
dvdt = squeeze(GetVar(statefile, diagfile, {'TOTVTEND', '(1)./86400'}, sliceD));

% ADVECTIVE TERMS
U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, sliceD);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, sliceD);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, sliceD);
Ux = Drv(500,U, 'x'); Uy = Drv(500,U, 'y'); Uz = Drv(metric, U, 'z');
Vx = Drv(500,V, 'x'); Vy = Drv(500,V, 'y'); Vz = Drv(metric,V, 'z');
UADV = U.*Ux + V.*Uy + W.*Uz;
VADV = U.*Vx + V.*Vy + W.*Vz;

% Coriolis Terms
Vm_Cori = GetVar(statefile, diagfile, {'Vm_Cori', '(1)'}, sliceD);
Um_Cori = GetVar(statefile, diagfile, {'Um_Cori', '(1)'}, sliceD);

% Calculate Pressure Gradients
dpdx = GetVar(statefile,diagfile, {'Um_dPHdx','(1)'},sliceD);
dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.81*(1)'},sliceEta));
dEtadx = Drv(500,dEtadx,'x');
dpdx = permute(repmat((dEtadx), [1, 1, 1, zL]), [1 2 4 3]) + dpdx;

dpdy = GetVar(statefile,diagfile, {'Vm_dPHdy','(1)'},sliceD);
dEtady = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.81*(1)'},sliceEta));  
dEtady = Drv(500,dEtady, 'y');
dpdy = permute(repmat((dEtady), [1, 1, 1, zL]), [1 2 4 3]) + dpdy;

UNC = dudt + UADV - Um_Cori -dpdx; %% XXX- Good place for a sanity check regarding closed budget...
VNC = dvdt + VADV - Vm_Cori -dpdy; %% XXX- Good place for a sanity check regarding closed budget...

%% BUOYANCY BUDGET
TtoB = 9.81.*2e-4;

dbdt = TtoB.*squeeze(GetVar(statefile, diagfile, {'TOTTTEND', '(1)./86400'}, sliceD));

b = TtoB.*squeeze(GetVar(statefile, diagfile, {'THETA', '(1)'}, sliceD));
bx = Drv(500, b, 'x'); by = Drv(500, b, 'y'); bz = Drv(metric, b, 'z');
BADV = U.*bx + V.*by + W.*bz;

BNC = dbdt + BADV;
% translate into cross-front coordinate.
theta = atan2(-squeeze(nanmean(nanmean(nanmean(by(:,:,1,:),4)))), -squeeze(nanmean(nanmean(nanmean(bx(:,:,1,:), 4))))) 
theta = atan2(-squeeze(by(:,:,:,:)), -squeeze(bx(:,:,:,:))) ;

ducdt = dvdt.*cos(pi/2-theta) + dudt.*cos(theta); % THETA calculated in MLIDiabaticPlot.m
UCADV = VADV.*cos(pi/2-theta) + UADV.*cos(theta);
UCCOR = Vm_Cori.*cos(pi/2-theta) + Um_Cori.*cos(theta);
UCPRES = dpdy.*cos(pi/2-theta) + dpdx.*cos(theta);
UCNC = VNC.*cos(pi/2-theta) + UNC.*cos(theta);


UACROSS = V.*cos(pi/2-theta) + U.*cos(theta);
UALONG = U.*cos(pi/2-theta) - V.*cos(theta);
DUAdAc = Drv(500, UALONG, 'x').*cos(theta) + Drv(500, UALONG, 'y').*cos(pi/2-theta);

CADV = UACROSS.*DUAdAc;
% DACROSS = Vy.*cos(pi/2-theta) + Uy.*cos(pi/2-theta) + Vx.*cos(theta)+Ux.*cos(theta);
% ADCROSS = UACROSS.*DACROSS;


DBACROSS = bx.*cos(theta) + by.*cos(pi/2-theta);  
EBF = squeeze(nanmean((UCCOR(ind,ind,:,:)+UCPRES(ind,ind,:,:)).*DBACROSS(ind,ind,:,:),4))./f0;
EBF = squeeze(nanmean(UCNC(ind,ind,:,:).*DBACROSS(ind,ind,:,:),4))./f0;
% EBFa = -squeeze(nanmean(nanmean(nanmean(UCNC.*DBACROSS, 4))))./f0;
EBFa = squeeze(nanmean(nanmean(nanmean((UCCOR+UCPRES).*DBACROSS, 4))))./f0;

% EBF = squeeze(nanmean(UACROSS(ind,ind,:,:).*DBACROSS(ind,ind,:,:),4))
Hkpp = squeeze(nanmean(GetVar(statefile,etanfile, {'KPPhbl', '(1)'}, sliceEta),4));
de = .4.*(9.81*2e-4*squeeze(25)./(1035*3994).*Hkpp).^(1/3)./f0;
% de = Hkpp;
EBF =  (-1+de(ind,ind)./Hkpp(ind,ind)).*EBF-1.2.*9.81*2e-4*squeeze(25)./((1035*3994.*Hkpp(ind,ind))) ;
EBFa = -squeeze(nanmean(nanmean((-1+de./Hkpp)))).*EBFa-squeeze(nanmean(nanmean(1.2.*9.81*2e-4*squeeze(25)./((1035*3994.*Hkpp))))) ;
% EBFa = -EBFa;
%%

figure
subplot(1,2,1)
plot(squeeze(nanmean(ducdt(ind,ind,:, :),4)), ztmp);
hold on
plot(squeeze(nanmean(UCADV(ind,ind,:, :),4)), ztmp);
plot(squeeze(nanmean(-UCCOR(ind,ind,:, :),4)), ztmp);
plot(squeeze(nanmean(-UCPRES(ind,ind,:, :),4)), ztmp);
plot(squeeze(nanmean(UCNC(ind,ind,:, :),4)), ztmp);
plot(squeeze(nanmean(-UCCOR(ind,ind,:, :)-UCPRES(ind,ind,:,:),4)), ztmp);
plot(squeeze(nanmean(CADV(ind,ind,:, :),4)), ztmp, '--');
% plot(squeeze(nanmean(CADV(ind,ind,:, :)-UCCOR(ind,ind,:,:)-UCPRES(ind,ind,:,:),4)), ztmp, '--');

% plot(squeeze(nanmean(DACROSS(ind,ind,:, :),4)), ztmp, 'r');

hold off

set(gca, 'ylim', [-75 -5]);
grid on
legend('dudt', 'UADV', 'fV', 'dpdx', 'NC')

subplot(1,2,2)
plot(squeeze(nanmean(dbdt(ind,ind,:,:),4)), ztmp);
hold on
plot(squeeze(nanmean(BADV(ind, ind,:,:),4)), ztmp);
plot(squeeze(nanmean(BNC(ind, ind,:,:),4)), ztmp)
% plot(squeeze(nanmean(dbdt(ind,ind,:,:)+BADV(ind,ind,:,:),4)), ztmp);
% plot(squeeze(nanmean(BNC(ind, ind,:,:),4))-EBF, ztmp)
plot(EBF, ztmp,'r')
hold off

set(gca, 'ylim', [-75 -5]);
grid on
legend('dbdt', 'BADV',  'NC')
%% AVERAGED
figure
subplot(1,2,1)
plot(squeeze(nanmean(nanmean(nanmean(ducdt(:,:,:, :),4)))), ztmp);
hold on
plot(squeeze(nanmean(nanmean(nanmean(UCADV(:,:,:, :),4)))), ztmp);
% plot(squeeze(nanmean(nanmean(nanmean(-UCCOR(:,:,:,:),4)))), ztmp);
% plot(squeeze(nanmean(nanmean(nanmean(-UCPRES(:,:,:, :),4)))), ztmp);
plot(squeeze(nanmean(nanmean(nanmean(UCNC(:,:,:, :),4)))), ztmp);
plot(squeeze(nanmean(nanmean(nanmean(-UCCOR(:,:,:, :)-UCPRES(:,:,:,:),4)))), ztmp, '--');
plot(squeeze(nanmean(nanmean(nanmean(CADV(:,:,:, :),4)))), ztmp);

set(gca, 'ylim', [-75 -5]);
grid on

hold off
subplot(1,2,2)
plot(squeeze(nanmean(nanmean(nanmean(dbdt(:,:,2:end, :),4)))), ztmp(2:end));
hold on
plot(squeeze(nanmean(nanmean(nanmean(BADV(:,:,:, :),4)))), ztmp);

plot(squeeze(nanmean(nanmean(nanmean(BNC(:,:,:, :),4)))), ztmp);
plot(EBFa(2), ztmp(2),'xr', 'LineWidth', 2, 'MarkerSize', 8)


% plot(squeeze(nanmean(DACROSS(ind,ind,:, :),4)), ztmp, 'r');

hold off
set(gca, 'ylim', [-75 ztmp(2)]);
grid on
%%
cl = [-1 1].*2e-6;

ti = floor(nt./2);
figure
subplot(2,2,1)
pcolor(1:nx, ztmp, squeeze(ducdt(ind,:,:,ti)).'); shading interp
colorbar
set(gca, 'clim', cl);

subplot(2,2,2)
pcolor(1:nx, ztmp, squeeze(UCADV(ind,:,:,ti)).'); shading interp
colorbar
set(gca, 'clim', cl);

subplot(2,2,3)
pcolor(1:nx, ztmp, squeeze(-UCCOR(ind,:,:,ti)-UCPRES(ind,:,:,ti)).'); shading interp
colorbar
set(gca, 'clim', cl);


subplot(2,2,4)
pcolor(1:nx, ztmp, squeeze(UCNC(ind,:,:,ti)).'); shading interp
colorbar
set(gca, 'clim', cl);
