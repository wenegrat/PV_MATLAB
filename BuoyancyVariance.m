bbar = nanmean(b);

bp  = b - repmat(bbar, [nx 1 1 1]);
vp  = V - repmat(nanmean(V), [nx 1 1 1]);
bbary = DPeriodic(bbar, dy, 'y');

wp = W - repmat(nanmean(W), [nx 1 1 1]);

%%
bbarz = Drv(metric(2,:,:,:), bbar, 'z');

vpbp = nanmean(vp.*bp);
wpbp = nanmean(wp.*bp);

EddyAdvectionY = vpbp.*bbary;
EddyAdvectionZ = wpbp.*bbarz;
EddyAdvection = EddyAdvectionY + EddyAdvectionZ;
EddyResidual = abs(bbary+1i.*bbarz).*vpbp;

Diffusion = nanmean(bp.*Q);%Note sign convention in Kuo et al. 2005

bvar = 1/2.*bp.*bp;
bvar = nanmean(bvar);

bvarx = nanmean(U).*DPeriodic(bvar, dx, 'x');
bvary = nanmean(V).*DPeriodic(bvar, dy, 'y');
gamma = squeeze(nanmean((bbary(:,:,2,:)))./nanmean((bbarz(:,:,2,:))));

[~, ~, ~, DdtBvar]= gradient(bvar, ts);
DdtBvar = - DdtBvar;


triplecorr = nanmean(1./2.*vp.*bp.*bp);
tcorry = DPeriodic(triplecorr, dy, 'y');
tcorrz = Drv(metric(2,:,:,:), nanmean(1./2.*wp.*bp.*bp), 'z');
tcorr = - tcorry-tcorrz;
%%

rhs = Diffusion + DdtBvar+tcorr+bvarx+bvary;
TotalDeriv = DdtBvar + bvarx + bvary;

yl =1:ny;

figure
plot(squeeze(nanmean(EddyAdvection(:,yl,2,:))));
hold on
plot(squeeze(nanmean(Diffusion(:,yl,2,:))));
plot(squeeze(nanmean(TotalDeriv(:,yl,2,:))));
plot(squeeze(nanmean(tcorr(:,yl,2,:))));
% plot(squeeze(nanmean(bvarx(:,:,2,:))));
% plot(squeeze(nanmean(bvary(:,:,2,:))));

plot(squeeze(nanmean(rhs(:,yl,2,:))), '--')
hold off

legend('Eddy Advection', 'Diffusion', 'DDt BVar', 'Triple Corr', 'SUM')

%%
plot(squeeze(nanmean(EddyAdvection(:,yl,2,:))));
hold on
plot(squeeze(nanmean(EddyResidual(:,yl,2,:))));
plot(squeeze(nanmean(EddyAdvection(:,yl,2,:) - EddyResidual(:,yl,2,:))));
hold off
%%
scatter(squeeze(nanmean(nanmean(outputFull.JBz(:,:,2,st:(st+nt-1))))), -squeeze(nanmean(nanmean(Q(:,:,2,:)))).*f0)
onetoone
cr = corr(squeeze(nanmean(nanmean(outputFull.JBz(:,:,2,st:(st+nt-1))))), -squeeze(nanmean(nanmean(Q(:,:,2,:)))).*f0);
title(num2str(cr))
