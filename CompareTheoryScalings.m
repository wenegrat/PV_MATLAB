%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare Theoretical Scalings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0 = 1e-4;
tslice = [1 599];
slice={0, 0, 0, tslice};
sliceEta={0,0,[1 1],tslice};

Q0 = ncread(etanfile, 'TFLUX');
THETA = GetVar(statefile, diagfile, {'THETA', '(1)'}, slice);
gradb = TtoB.*(DPeriodic(THETA, dx, 'x').^2 + DPeriodic(THETA, dy, 'y').^2).^(1/2);
Hbl = GetVar(statefile, etanfile, {'KPPhbl', '(1)'}, sliceEta);

Zfull = permute(repmat(Z, [1,nx, ny, nt]), [2 3 1 4]);
[nx, ny, nz, nt] = size(gradb);
mask = Zfull>-repmat(Hbl, [1, 1, nz, 1]);

masknan = double(mask);
masknan(masknan<1) = NaN;
gradba = squeeze(nanmean(nanmean(nanmean(gradb.*masknan))));
gradban = gradba./gradba(1); % ONLY APPROXIMATE
gradba = sqrt(squeeze(nanmean(nanmean(nanmean((gradb.^2).*masknan)))));
Hbla = squeeze(nanmean(nanmean(Hbl)));
B0 = 9.81*2e-4*squeeze(Q0(1,1,1,1))./(1035*3994);


gammaFD = gradba.^2.*Hbla.^(4/3)./(f0^2.*abs(B0).^(2/3));
% gammaFD = gradba.^2./(f0^2.*abs(B0).^(2/3));
% vargradb = squeeze(nanvar(nanvar(nanvar(gradb.*masknan))));
%%
%AreaIntegrateJTerms
slice={0, 0, [1 5], tslice};
sliceEta={0,0,[1 1],tslice};
nz = 5;
CalculateQTerms;
isoT = [0 100]; %Needed for plot below
mask = ones(nx, ny, nz, tslice(end)-tslice(1)+1);
vol = gridvol.*squeeze(sum(sum(sum(mask))));
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);

[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
%%
gammaDir = abs(dJFdt./dJBdt);

%%
if ~exist('fgammaFD');
fgammaFD = figure;
else
    figure(fgammaFD);
end
hold on
plot(time, gammaFD/100, 'LineWidth',2);
plot(time, gammaDir, 'LineWidth',2 ,'LineStyle', '--');
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

