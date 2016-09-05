%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare Theoretical Scalings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0 = 1e-4;
tslice = [1 599];
slice={0, 0, 0, tslice};
sliceEta={0,0,[1 1],tslice};

% Q0 = ncread(etanfile, 'TFLUX', slice);
Q0 = GetVar(statefile, etanfile, {'TFLUX', '(1)'}, sliceEta);
THETA = GetVar(statefile, diagfile, {'THETA', '(1)'}, slice);
gradb = TtoB.*(DPeriodic(THETA, dx, 'x').^2 + DPeriodic(THETA, dy, 'y').^2).^(1/2);
ZETA = GetVar(statefile, extrafile, {'momVort3', '(1)'}, sliceEta);
Hbl = GetVar(statefile, etanfile, {'KPPhbl', '(1)'}, sliceEta);
ghat = GetVar(statefile, 'kppdiags.nc', {'KPPg_TH',['(1)/',num2str(dx*dy)]}, slice);
bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
magbgrad = bx.^2 + by.^2;
K = GetVar(statefile, 'kppdiags.nc', {'KPPdiffT', '(1)'}, slice);
Zfull = permute(repmat(Z, [1,nx, ny, nt]), [2 3 1 4]);
[nx, ny, nz, nt] = size(gradb);
mask = Zfull>-repmat(Hbl, [1, 1, nz, 1]);

masknan = double(mask);
masknan(masknan<1) = NaN;
gradba = squeeze(nanmean(nanmean(nanmean(gradb.*masknan))));
gradban = gradba./gradba(1); % ONLY APPROXIMATE
gradba = sqrt(squeeze(nanmean(nanmean(nanmean((gradb.^2).*masknan)))));
gradbvar = (squeeze(nanvar(nanvar(nanvar((gradb.*masknan))))));
%%
Hbla = squeeze(nanmean(nanmean(Hbl)));
B0 = 9.81*2e-4*squeeze(Q0(1,1,1,1))./(1035*3994);

% rho = 1035*(1-2e-4.*(THETA(:,:,1,:)-16.5));
% B0 = 9.81*2e-4*squeeze(Q0(1,1,1,1))./(rho*3994);

% JfT = Hbla.^(4/3).*abs(B0).^(1/3).*gradba.*2;
% JfT = Hbla.^(1/3).*abs(B0).^(1/3).*gradbvar;

Jftot = squeeze(gradb(:,:,1,:)).*squeeze(gradb(:,:,1,:)).*squeeze(Hbl).^(4/3).*abs(squeeze(B0)).^(1/3);
JBtot = B0.*(f0+ZETA)./Hbl;
Jftota = squeeze(nanmean(nanmean(Jftot))).*vol;
Jbtota = squeeze(nanmean(nanmean(JBtot)));
%%
JBtotp = B0.*f0./(Hbl); %Corr = 0.73
% JBtotp = (f0+ZETA).*(B0-K(:,:,2,:).*bz(:,:,2,:)+TtoB.*ghat(:,:,2,:)); % All 3 terms corr = 0.95
% JBtotp = (f0+ZETA).*(B0-K(:,:,2,:).*bz(:,:,2,:)); %Only Surf and Local corr = 0.89
JBtotp = (f0+ZETA).*(B0+TtoB.*ghat(:,:,2,:)); % Only surf and nonlocal corr=0.68,

JBtotp = f0.*(B0 - 0.06.*magbgrad(:,:,2,:).*Hbl.^2./f0)./Hbl; %Surf and <w'b'> Fox-Kemper corr = 0.92.
JBtotp = f0.*(magbgrad(:,:,2,:).*Hbl.^2./f0)./Hbl; % Only <w'b'> corr = 0.95
% JBtotp = (f0+ZETA).*(B0-Hbl.*abs(B0.*Hbl).^(1/3).*bz(:,:,2,:)+TtoB.*ghat(:,:,2,:));

% JBtotp = K(:,:,2,:).*bz(:,:,2,:);
% JBtotp = -(f0).*(squeeze((abs(B0).*Hbl).^(1/3)).*TtoB.*squeeze(THETA(:,:,1,:))-B0)./squeeze(Hbl);
% JBtotp = -f0*(squeeze((abs(B0).*Hbl).^(1/3)).*TtoB.*squeeze(THETA(:,:,1,:)))./squeeze(Hbl);
% JBtotp = f0*B0.*(1-Hbl)./Hbl;
% JBtotp = abs(f0*(Q0 - Hbl.*(Hbl.*squeeze(B0)).^(1/3).*(bz(:,:,2,:))));
Jbtotpa = squeeze(nanmean(nanmean(JBtotp))).*vol;
% Jbtotpa = f0*B0./3;
% gammaFD = gradba.^2.*Hbla.^(4/3)./(f0^2.*abs(B0).^(2/3));
% 
% % gammaFD = gradba.^2.*Hbla.^(5/3)./(f0.^(3/2).*abs(B0).^(5/6));
% gammaFD = gradbvar.*Hbla.^(4/3)./(f0^2.*abs(B0).^(2/3));
gammaFD = Jftota./Jbtota;
% gammaFD = gradba.^2./(f0^2.*abs(B0).^(2/3));
% vargradb = squeeze(nanvar(nanvar(nanvar(gradb.*masknan))));
%%
% Plot the 3 components of Diabatic term
figure
plot(squeeze(nanmean(nanmean(B0.*f0.*ones(size(Hbl))))));
hold on
plot(squeeze(nanmean(nanmean(f0.*(K(:,:,2,:).*bz(:,:,2,:))))));
plot(squeeze(nanmean(nanmean(f0.*TtoB.*ghat(:,:,2,:)))));
% plot(squeeze(nanmean(nanmean(f0.*(K(:,:,2,:).*bz(:,:,2,:))+ f0.*TtoB.*ghat(:,:,2,:)))),'--');

hold off
% Note that ghat and K term are highly correlated, both dependent on HBL
%%
scalen2 = squeeze(nanmean(nanmean(magbgrad(:,:,2,:).*Hbl.^2./f0)));
plotyy(1:599, squeeze(nanmean(nanmean(f0.*(K(:,:,2,:).*bz(:,:,2,:))))),1:599, f0.*0.06.*scalen2)

%%
plotyy(1:599, squeeze(nanmean(nanmean(K(:,:,2,:)))), 1:599, squeeze(nanmean(nanmean(Hbl.*abs(B0.*Hbl).^(1/3).*3))));
title(num2str(corr(squeeze(nanmean(nanmean(K(:,:,2,:)))), squeeze(nanmean(nanmean(Hbl.*abs(B0.*Hbl).^(1/3).*3))))))
%%
% %AreaIntegrateJTerms
slice={0, 0, [1 5], tslice};
sliceEta={0,0,[1 1],tslice};
nz = 5;
% CalculateQTerms;
isoT = [0 100]; %Needed for plot below
mask = ones(nx, ny, nz, tslice(end)-tslice(1)+1);
vol = gridvol.*squeeze(sum(sum(sum(mask))));
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);

[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);

gammaDir = abs(dJFdt./dJBdt);
%%
subplot(3,1,1)
% indj = 30; indk = 76;
% plotyy(time, abs(squeeze(JFz(indj, indk, 2,:))), time, squeeze(JfTtot(indj, indk,:)));
[ax, h1, h2] = plotyy(time, -dJFdt, time, Jftota);
set(h1, 'LineWidth', 2);
set(h2, 'LineWidth', 2);
grid on
title(num2str(corr(dJFdt, Jftota)))

subplot(3,1,2)
[ax, h1, h2]=plotyy(time, -dJBdt, time, -Jbtotpa);
set(h1, 'LineWidth', 2);
set(h2, 'LineWidth', 2);
grid on
hold on
% plot(time, Jbtotpa, '--');
hold off
title(num2str(corr(dJBdt, Jbtotpa)))

subplot(3,1,3)
% plotyy(time, abs(dJFdt./dJBdt), time, abs(Jftota./Jbtota));
plot(time, abs(dJFdt./dJBdt), 'LineWidth',2);
hold on
plot(time, abs(Jftota./Jbtotpa), 'LineWidth', 2);
hold off
grid on
title(num2str(corr(abs(dJFdt./dJBdt), abs(Jftota./Jbtotpa))))
set(gca, 'ylim', [0 10]);
%%
if ~exist('fgammaFD');
fgammaFD = figure;
else
    figure(fgammaFD);
end
hold on
plot(time, gammaFD, 'LineWidth',2);
% plot(time, gammaDir, 'LineWidth',2 ,'LineStyle', '--');
hold off
title('\gamma_{FD}');

grid on
xlabel('Days');
set(gcf, 'Color','w');

%%
if ~exist('fkpphbl');
fkpphbl= figure;
else
    figure(fkpphbl);
end
hold on
plot(time, Hbla, 'LineWidth',2);
hold off
title('h_{BL}');

grid on
xlabel('Days');
set(gcf, 'Color','w');

%%
if ~exist('fhgb');
fhgb = figure;
else
    figure(fhgb);
end
hold on
plot(time, gradba, 'LineWidth',2);
hold off
title('\nabla b');

grid on
xlabel('Days');
set(gcf, 'Color','w');

%%
if ~exist('fhgvar');
fhgvar = figure;
else
    figure(fhgvar);
end
hold on
plot(time, gradbvar, 'LineWidth',2);
hold off
title('Var b');

grid on
xlabel('Days');
set(gcf, 'Color','w');