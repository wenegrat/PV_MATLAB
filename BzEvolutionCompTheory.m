%% Rate of Change of Stratification
statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc'; extrafile = 'extra.nc';
kppfile = 'kppdiags.nc';


hl = 36;
Z = ncread(statefile, 'Z');
Zl = ncread(statefile, 'Zl');
nx = 150; ny =200; 


dx = 500;
dy = 500;
ts = 7200;
st = 50;
zl = 45;
nt = 100;
Z = Z(1:zl);
nz = zl;
metric = permute(repmat(Z,  [1, nx, ny, nt]), [2 3 1 4]);

slice = {0, 0, [1 zl], [st st+nt-1]};
Ho = abs(Z(hl));
b = GetVar(statefile, diagfile, {'b', '(1)'},slice);

%% LHS
% Note Assuming constant H...
bzbar = squeeze(nanmean(nanmean((b(:,:,1,:) - b(:,:,hl,:))./(Ho))));
nzt = gradient(bzbar, ts);

%% RHS
% Fox-Kemper Parameterization
Ce = 0.06;

hkpp = GetVar(statefile, etanfile, {'KPPhbl', '(1)'}, {slice{1},slice{2}, [1 1], slice{4}});

bax = nanmean(b);
baxz = nanmean(bax(:,:,1:hl,:), 3);
bx = DPeriodic(b, dx, 'x');
by = DPeriodic(b, dy, 'y');
magb = abs(bx+1i*by);
%%
% Average in x,y, and over MLD
% ms = double(metric > -repmat(hkpp, [1 1 nz 1]));
% ms(ms==0)  = NaN;
bxbar = squeeze(nanmean(nanmean(nanmean(abs(bx(:,:,1:hl,:))))));
bybar = squeeze(nanmean(nanmean(nanmean(abs(by(:,:,1:hl,:))))));

% bybar = squeeze(nanmean(nanmean(nanmean(abs(by).^2))));
magba = squeeze(nanmean(nanmean(nanmean(magb(:,:,1:hl,:)))));

zt = repmat(Z, [1 nt]);
hz = repmat(squeeze(nanmean(nanmean(hkpp))), [1 nz]).';
dmudz = - 8*(hz + 2*zt).*(13*hz.^2 +20 * hz.*zt + 20.*zt.^2)./(21*hz.^4);
ms = zt > - hz;
dmudz = ms.*dmudz;


dmudzA = (dmudz(1,:) - dmudz(hl,:))./Ho;
% dmudzA = 1;
EddyFlux = -Ce.*squeeze(nanmean(nanmean(hkpp)).^2).*(magba.^2)./f0.*dmudzA.';

b0 = 9.81*2e-4*squeeze(25)./(1035*3994);

Dsurf = -b0./squeeze(nanmean(nanmean(hkpp)))./Ho;
Deddy = 0.05.*squeeze(nanmean(nanmean(hkpp))).*(bxbar.^2 + bybar.^2)./f0./Ho;
%%
figure
plot(1:nt, nzt);
hold on
plot(1:nt, EddyFlux);
plot(1:nt, Deddy);
plot(1:nt, Dsurf);
plot(1:nt, EddyFlux+Dsurf+Deddy, '--')
hold off

legend('dNzdt','EddyFlux','DEddy', 'DSurf');

%%
% scatter(nzt, EddyFlux)












