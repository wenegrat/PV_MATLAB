%% Extra Diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) EKE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Fox-Kemper Buoyancy scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ce = 0.06;
f0 = 1e-4;

Q0 = ncread(etanfile, 'TFLUX');
THETA = GetVar(statefile, diagfile, {'THETA', '(1)'}, slice);
gradb = TtoB.*(DPeriodic(THETA, dx, 'x').^2 + DPeriodic(THETA, dy, 'y').^2).^(1/2);

[nx, ny, nz, nt] = size(gradb);
deltaCrit = 0.01./(-2e-4*1035); % Convert to a delta T criteria.
tprime = THETA - repmat(THETA(:,:,1,:), [1, 1, nz, 1]);
mask = tprime>deltaCrit;
Zfull = permute(repmat(Z, [1,nx, ny, nt]), [2 3 1 4]);
mld = squeeze(min(Zfull.*mask, [], 3));

%%
masknan = double(mask);
masknan(masknan<1) = NaN;
gradba = squeeze(nanmean(nanmean(nanmean(gradb.*masknan))));
vargradb = squeeze(nanvar(nanvar(nanvar(gradb.*masknan))));
mlda = squeeze(nanmean(nanmean(mld)));

psib = ce.*abs(gradba).^2.*mlda.^2./f0;

B0 = 9.81*2e-4*squeeze(Q0(1,1,1,1))./(1035*3994);

%%
% Make Buoyancy Flux Plot
if ~exist('fhpsi');
fhpsi = figure;
else
    figure(fhpsi);
end
hold on
plot(time, psib./B0, 'LineWidth',2);
hold off
title('<w''b''>_{FK}/B_0');

grid on
xlabel('Days');
set(gcf, 'Color','w');
%%
% Make MLD  PLOT
if ~exist('fhmld');
fhmld = figure;
else
    figure(fhmld);
end
hold on
plot(time, mlda, 'LineWidth',2);
hold off
title('MLD');

grid on
xlabel('Days');
set(gcf, 'Color','w');

%%
% Make Buoyancy Flux Plot
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