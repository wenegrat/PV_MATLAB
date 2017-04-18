%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check subinertial mom budget for dominate balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc'; extrafile = 'extra.nc';
kppfile = 'kppdiags.nc';
nt = 48;
st = 48;
% st = 120;
tslice = [st st+nt-1];
slice={0, 0, 0, tslice};

nx = 150; ny = 200; nz = 50;
nx = 3;
sizes = [nx, ny, nz, nt];
dx = 500; dy = 500; dz = 0;
[DPDY, VADV, VCORI, VTEND, IFRICV, TY, Vx, Vy, Vz , TYY]= returnMomTerms(diagfile, statefile, etanfile,extrafile, sizes, slice, dx, dy,dz );
Z = ncread(statefile, 'Z');
%%

dpdy = nanmean(DPDY, 4);
vadv = nanmean(VADV, 4);
vcori = nanmean(VCORI, 4);
vtend = nanmean(VTEND, 4);
ifricv = nanmean(IFRICV, 4);
ty = nanmean(TY, 4);
% ty = abs(ty); %xxx-think about this;
lhs =  vtend - dpdy - vadv - vcori;
ageo = -(vcori + dpdy);
vx = nanmean(Vx, 4);
vy = nanmean(Vy, 4);
vz = nanmean(Vz, 4);

dpdya = nanmean(nanmean(dpdy, 2));
vadva = nanmean(nanmean(vadv, 2));
vcoria = nanmean(nanmean(vcori,2 ));
vtenda = nanmean(nanmean(vtend,2));
ifrica = nanmean(nanmean(ifricv, 2));
ageoa = nanmean(nanmean(ageo, 2));
vxa = nanmean(nanmean(vx, 2));
vya = nanmean(nanmean(vy, 2));
vza = nanmean(nanmean(vz, 2));

dpdyaw = nanmean(nanmean(dpdy.*ty));
vadvaw = nanmean(nanmean(vadv.*ty));
vcoriaw = nanmean(nanmean(vcori.*ty));
vtendaw = nanmean(nanmean(vtend.*ty));
ifricaw = nanmean(nanmean(ifricv.*ty));
ageoaw = nanmean(nanmean(ageo.*ty));
%%2
% Figures

xi = 2; yi = 2;

figure
subplot(1,3,1);
plot(squeeze(ifricv(xi, yi, :,:)), Z);
hold on
plot(squeeze(vtend(xi, yi,:,:)), Z, '-');
plot(squeeze(-vadv(xi, yi,:,:)), Z, '-');
% plot(squeeze(-vcori(xi, yi,:,:)), Z, '-');
% plot(squeeze(-dpdy(xi, yi,:,:)), Z, '-');
plot(squeeze(ageo(xi,yi,:,:)), Z, 'LineWidth', 2, 'LineStyle', '--');
plot(squeeze(vtend(xi,yi,:,:) +ageo(xi,yi,:,:) -vadv(xi,yi,:,:) ), Z, 'r--');
hold off

legend('FRIC', 'TEND', 'ADV', 'AGEO')

% plot(-squeeze(-vtendaw + vadvaw +ifricaw), Z); %Ageostrophic advection of Ty
% hold on
% plot(-squeeze(-vtendaw), Z); 
% plot(-squeeze(vadvaw), Z); 
% plot(-squeeze(ifricaw), Z); 
% hold off
% legend('vTy', 'VTEND', 'VADV', 'FRIC');

subplot(1,3,2);
plot(squeeze(ifrica), Z);
hold on
plot(squeeze(vtenda), Z, '-');
plot(squeeze(-vadva), Z, '-');
% plot(squeeze(-vcori(xi, yi,:,:)), Z, '-');
% plot(squeeze(-dpdy(xi, yi,:,:)), Z, '-');
plot(squeeze(ageoa), Z, 'LineWidth', 2, 'LineStyle', '--');
plot(squeeze(ageoa - vadva), Z, ':')
% plot(squeeze(vtenda - vadva + ageoa - ifrica), Z);
% plot(squeeze(-vxa), Z, '-');
% plot(squeeze(-vya), Z, '-');
% plot(squeeze(-vza), Z, '-');

hold off

legend('FRIC', 'TEND', 'ADV', 'AGEO')

subplot(1,3,3);
plot(squeeze(ifricaw), Z);
hold on
plot(squeeze(vtendaw), Z, '-');
plot(squeeze(-vadvaw), Z, '-');
% plot(squeeze(-vcori(xi, yi,:,:)), Z, '-');
% plot(squeeze(-dpdy(xi, yi,:,:)), Z, '-');
plot(squeeze(ageoaw), Z, 'LineWidth', 2, 'LineStyle', '--');
plot(squeeze(ageoaw - vadvaw), Z, ':');
hold off

legend('FRIC', 'TEND', 'ADV', 'AGEO')