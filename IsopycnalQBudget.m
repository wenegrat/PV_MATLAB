%% 
tic
% GENERATE ISOPYCNAL BUDGET
  clc; clear all; %close all;
 statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
% Parameters            
% dx = 500; dy = dx; dz = 3;
% nx = 48; ny=48; nz=200;
% ts = 3600;

dx = 250; dy = dx; dz = 3;
nx = 96; ny=144; nz=100;
ts = 3600;
% 
% dx = 250; dy = dx; dz = 3;
% nx = 96; ny=192; nz=200;
% ts = 3600;

TtoB = 9.81.*2e-4;
tslice = [10 200];
tslice = [230 309];
 tslice = [1 719];
slice={0, 0, 0, tslice};%100 120
sliceEta={0,0,[1 1],tslice};%251 271


%%
% CALCULATE TERMS
Q = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);
Qdir = Q;
JBx = Q;
JBy = Q;
JBz = Q;
JFx = Q;
JFz = Q;
JFy = Q;
JAx = Q;
JAy = Q;
JAz = Q;
JFzN = Q;
JBzN = Q;
mask = Q;
FricDiv = Q;
AdvDiv = Q;
DiaDiv = Q;
THETA  = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);
resid = Q;
uterm = Q;
bterm = Q;
% remaining = (tslice(end)-tslice(1)+1);
incc =1;
inc = incc-1;

ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, nx, ny, incc]), [2 3 1 4]);
% for i=1:incc:(tslice(end)-tslice(1)+1-inc)
%     disp(num2str(i));
% %     remaining = @(x)    remaining -1
% %     disp(['Remaining: ', num2str(remaining), '    (Processing: ', num2str(i),')'])
%    slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1+inc]};
%    [Q(:,:,:,i:i+inc), Qdir(:,:,:,i:i+inc), JAx(:,:,:,i:i+inc), JAy(:,:,:,i:i+inc), JAz(:,:,:,i:i+inc), ...
%        JFx(:,:,:,i:i+inc), JFy(:,:,:,i:i+inc), JFz(:,:,:,i:i+inc), JBx(:,:,:,i:i+inc), JBy(:,:,:,i:i+inc), JBz(:,:,:,i:i+inc), JFzN(:,:,:,i:i+inc), JBzN(:,:,:,i:i+inc)] ...
%        = calcQBudgetD(diagfile, statefile, etanfile, [nx, ny, incc], slicetemp, dx, dy);
% 
%    FricDiv(:,:,:,i:i+inc) = Drv(dx, JFx(:,:,:,i:i+inc), 'x') + Drv(dy, JFy(:,:,:,i:i+inc), 'y') + Drv(metric, JFz(:,:,:,i:i+inc), 'z');
% %    AdvDiv(:,:,:,i:i+inc) = Drv(dx, JAx(:,:,:,i:i+inc), 'x') + Drv(dy, JAy(:,:,:,i:i+inc), 'y') + Drv(metric, JAz(:,:,:,i:i+inc), 'z');
%    AdvDiv(:,:,:,i:i+inc) = DPeriodic(JAx(:,:,:,i:i+inc), dx, 'x') + DPeriodic(JAy(:,:,:,i:i+inc),dy, 'y') + Drv(metric, JAz(:,:,:,i:i+inc), 'z');
% 
%    DiaDiv(:,:,:,i:i+inc) = Drv(dx, JBx(:,:,:,i:i+inc), 'x') + Drv(dy, JBy(:,:,:,i:i+inc), 'y') + Drv(metric, JBz(:,:,:,i:i+inc), 'z');
%    THETA(:,:,:,i:i+inc) = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);
% %    [resid(:,:,:,i:i+inc), uterm(:,:,:,i:i+inc), bterm(:,:,:,i:i+inc)] = AssessQCancellation(statefile, diagfile, slicetemp);
% end

%Do end bit
% i = i +incc;
% slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(2)]};
% metric = permute(repmat(ztmp, [1, nx, ny, slicetemp{4}(2)-slicetemp{4}(1)+1]), [2 3 1 4]);
% 
%    [Q(:,:,:,i:end), Qdir(:,:,:,i:end), JAx(:,:,:,i:end), JAy(:,:,:,i:end), JAz(:,:,:,i:end), ...
%        JFx(:,:,:,i:end), JFy(:,:,:,i:end), JFz(:,:,:,i:end), JBx(:,:,:,i:end), JBy(:,:,:,i:end), JBz(:,:,:,i:end), JFzN(:,:,:,i:end), JBzN(:,:,:,i:end)] ...
%        = calcQBudgetD(diagfile, statefile, etanfile, [nx, ny, slicetemp{4}(2)-slicetemp{4}(1)+1], slicetemp, dx, dy);
% 
%    FricDiv(:,:,:,i:end) = Drv(dx, JFx(:,:,:,i:end), 'x') + Drv(dy, JFy(:,:,:,i:end), 'y') + Drv(metric, JFz(:,:,:,i:end), 'z');
% %    AdvDiv(:,:,:,i:i+inc) = Drv(dx, JAx(:,:,:,i:i+inc), 'x') + Drv(dy, JAy(:,:,:,i:i+inc), 'y') + Drv(metric, JAz(:,:,:,i:i+inc), 'z');
%    AdvDiv(:,:,:,i:end) = DPeriodic(JAx(:,:,:,i:end), dx, 'x') + DPeriodic(JAy(:,:,:,i:end),dy, 'y') + Drv(metric, JAz(:,:,:,i:end), 'z');
% 
%    DiaDiv(:,:,:,i:end) = Drv(dx, JBx(:,:,:,i:end), 'x') + Drv(dy, JBy(:,:,:,i:end), 'y') + Drv(metric, JBz(:,:,:,i:end), 'z');
%    THETA(:,:,:,i:end) = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);

   
parfor i=1:1:(tslice(end)-tslice(1)+1-inc)
    disp(num2str(i));
%     remaining = @(x)    remaining -1
%     disp(['Remaining: ', num2str(remaining), '    (Processing: ', num2str(i),')'])
   slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1+inc]};
   [Q(:,:,:,i), Qdir(:,:,:,i), JAx(:,:,:,i), JAy(:,:,:,i), JAz(:,:,:,i), ...
       JFx(:,:,:,i), JFy(:,:,:,i), JFz(:,:,:,i), JBx(:,:,:,i), JBy(:,:,:,i), JBz(:,:,:,i), JFzN(:,:,:,i), JBzN(:,:,:,i)] ...
       = calcQBudgetD(diagfile, statefile, etanfile, [nx, ny, incc], slicetemp, dx, dy);

   FricDiv(:,:,:,i) = Drv(dx, JFx(:,:,:,i), 'x') + Drv(dy, JFy(:,:,:,i), 'y') + Drv(metric, JFz(:,:,:,i), 'z');
%    AdvDiv(:,:,:,i:i+inc) = Drv(dx, JAx(:,:,:,i:i+inc), 'x') + Drv(dy, JAy(:,:,:,i:i+inc), 'y') + Drv(metric, JAz(:,:,:,i:i+inc), 'z');
   AdvDiv(:,:,:,i) = DPeriodic(JAx(:,:,:,i), dx, 'x') + DPeriodic(JAy(:,:,:,i),dy, 'y') + Drv(metric, JAz(:,:,:,i), 'z');

   DiaDiv(:,:,:,i) = Drv(dx, JBx(:,:,:,i), 'x') + Drv(dy, JBy(:,:,:,i), 'y') + Drv(metric, JBz(:,:,:,i), 'z');
   THETA(:,:,:,i) = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);
%    [resid(:,:,:,i:i+inc), uterm(:,:,:,i:i+inc), bterm(:,:,:,i:i+inc)] = AssessQCancellation(statefile, diagfile, slicetemp);
end


B0 = ncread(etanfile, 'TFLUX');
X = ncread(statefile, 'X');
Y = ncread(statefile, 'Y');
Z = ncread(statefile, 'Z');
T = ncread(statefile, 'T');

%%
% Generate Mask
mask = zeros(nx, ny, nz, tslice(end)-tslice(1)+1);
isoT = [16.55 16.9];
isoT = [ 17 17.25];
isoT = [16.5 16.75];
isoT = [16.25 16.5];
 isoT = [1 100];
% isoT = [16.1 16.3];
% isoT = [3 16.55];
[nx, ny, nz, nt] = size(THETA);
for i=1:nt;
%    slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};
   mask(:,:,:,i) = (THETA(:,:,:,i)>isoT(1)) & (THETA(:,:,:,i)<isoT(2));
end 

% mask(:,end-2:end,:,:) = 0; %wall effects that might be in Q
% mask(:,:,1:2,:) = 0; % Ad hoc removal of part of Q to compensate for calc of J vectors.
%%
tind = floor((tslice(2)-tslice(1))/2);
% 
% [c, h] = contourf(squeeze(THETA(:,:,1,tind)).', linspace(15, 20, 100)); shading interp
% set(h, 'edgecolor','none')
% xlabel('x'); ylabel('y');
% hold on
% contour(squeeze(THETA(:,:,1,tind)).', isoT, 'k')
% set(gca, 'clim', [16 17])
% colorbar
%%
% Mask by depth
%zl = [1 20];
%for i=1:(tslice(end)-tslice(1)+1)
%   mask(:,:,:,i) = ones(size(THETA(:,:,:,i)));
%   mask(:,:, 1:zl(1),i) = 0;
%   mask(:,:, zl(2):end, i) = 0;
%end
%%
% Calculate divergences
% Integrate in Volume.
disp('Integrate Over Volume');

%Params
[~, ~, ~, nt] = size(Q);
t = (1:nt).*ts;
gridvol = dx.*dy.*dz;
vol = dx.*dy.*dz.*squeeze(sum(sum(sum(mask))));
% vol =1;

% Calculate Fluxes as integrals over divergences.
%Friction
Frici = FricDiv;%.*mask;
Frici(~isfinite(Frici)) = 0;
Advi = AdvDiv;%.*mask;
Advi(~isfinite(Advi)) = 0;
Diai = DiaDiv;%.*mask;
Diai(~isfinite(Diai))=0;
Frici = cumtrapz(t, Frici, 4);
Frici = Frici.*mask;
% Fric = squeeze(nansum(nansum(nansum(Frici)))).*gridvol;
Fric = squeeze(trapz(trapz(trapz(Frici)))).*gridvol;
% Fric = cumtrapz(t, Fric);
Frict = Fric./vol;
Advi = cumtrapz(t, Advi, 4);
Advi = Advi.*mask;
Adv = squeeze(nansum(nansum(nansum(Advi)))).*gridvol;
% Adv = cumtrapz(t,Adv);
Advt = Adv./vol;
Diai = cumtrapz(t, Diai, 4);
Diai = Diai.*mask;
Dia = squeeze(nansum(nansum(nansum(Diai)))).*gridvol;
Fric = squeeze(trapz(trapz(trapz(Frici)))).*gridvol;

% Dia = cumtrapz(t,Dia);
Diat = Dia./vol;

%Calculate as surface fluxes
zlt = 2;
zl =1;
yln = 0;
yls = 1;
withsides = false;
%Original Flux Calc
Frics = squeeze(nansum(nansum(JFz(:,yls:end-yln,zlt,:).*mask(:,yls:end-yln,zlt,:))) - nansum(nansum(JFz(:,yls:end-yln,end-zl+1,:).*mask(:,yls:end-yln,end-zl+1,:)))).*dx.*dy;
% Frics = squeeze(nansum(nansum(JFz(:,yls:end-yln,zlt,:).*mask(:,yls:end-yln,zlt,:)))).*dx.*dy;

if withsides; Frics = Frics+squeeze(nansum(nansum(JFy(:,end-yln,:,:).*mask(:,end-yln,:,:))) - nansum(nansum(JFy(:,yls,:,:).*mask(:,yls,:,:)))).*dz.*dx; end;
Fricst = cumtrapz(Frics).*ts./vol;
% Fricst = Fricst - Fricst(1);

Dias = squeeze(nansum(nansum(JBz(:,yls:end-yln,zlt,:).*mask(:,yls:end-yln,zlt,:))) - nansum(nansum(JBz(:,yls:end-yln,end-zl+1,:).*mask(:,yls:end-yln,end-zl+1,:)))).*dx.*dy;
% Dias = squeeze(nansum(nansum(JBz(:,yls:end-yln,zlt,:).*mask(:,yls:end-yln,zlt,:))) ).*dx.*dy;

if withsides; Dias = Dias+squeeze(nansum(nansum(JBy(:,end-yln,:,:).*mask(:,end-yln,:,:))) - nansum(nansum(JBy(:,yls,:,:).*mask(:,yls,:,:)))).*dz.*dx; end;
Diast = cumtrapz(Dias).*ts./vol;
% Diast = Diast - Diast(1);

FricsN = squeeze(nansum(nansum(JFzN(:,yls:end-yln,zlt,:).*mask(:,yls:end-yln,zlt,:))) - nansum(nansum(JFzN(:,yls:end-yln,end-zl+1,:).*mask(:,yls:end-yln,end-zl+1,:)))).*dx.*dy;
FricstN = cumtrapz(FricsN).*ts./vol;
DiasN = squeeze(nansum(nansum(JBzN(:,yls:end-yln,zlt,:).*mask(:,yls:end-yln,zlt,:))) - nansum(nansum(JBzN(:,yls:end-yln,end-zl+1,:).*mask(:,yls:end-yln,end-zl+1,:)))).*dx.*dy;
DiastN = cumtrapz(DiasN).*ts./vol;

%Alternate Flux Calc
% JINT = cumtrapz(t, JFz, 4);
% % Frics = squeeze(nansum(nansum(JINT(:,yl:end-yl, zl,:).*mask(:,yl:end-yl,zl,:))) - nansum(nansum(JINT(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dy.*dy;
% Frics = squeeze(nansum(nansum(JINT(:,yls:end-yln, zlt,:).*mask(:,yls:end-yln,zlt,:))) ).*dy.*dy;
% 
% Fricst = Frics./vol;
% 
% BINT = cumtrapz(t, JBz, 4);
% % Dias = squeeze(nansum(nansum(BINT(:,yl:end-yl, zl,:).*mask(:,yl:end-yl,zl,:))) - nansum(nansum(BINT(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dy.*dy;
% Dias = squeeze(nansum(nansum(BINT(:,yls:end-yln, zlt,:).*mask(:,yls:end-yln,zlt,:))) ).*dy.*dy;
% 
% Diast = Dias./vol;


%Integrate dQdt from Mom equations (should follow same order as J vectors).
Qi = Qdir.*mask;
% Qi = Qdir;
Qi(~isfinite(Qi)) = 0;
Qi = cumtrapz(t, Qi, 4);
Qta = squeeze(nansum(nansum(nansum(Qi)))).*gridvol; %This is volume integral of dQ/dt
Qda = Qta./vol; %Time Integral of Vol Integral of dQ/dt.
Qda = Qda - Qda(1);


% Calculate from Q term.
Qi = Q;
Qi(~isfinite(Qi)) = 0;
Qi = Qi.*mask;
% Qi = Qi - repmat(Qi(:,:,:,1), [1 1 1 nt]);
Qa = squeeze(nansum(nansum(nansum(Qi)))).*gridvol;
Qt = gradient(Qa, ts);
Qa = Qa-Qa(1);
 Qa = Qa./vol;%

% Qa = Qa./1035
%%
% residual = cumtrapz(t, uterm+bterm, 4);
% residual = squeeze(nansum(nansum(nansum(residual.*mask)))).*gridvol;
% residual = residual./vol;
% residavg = residual;
%%
% Make Time Series Figure of DeltaQ and Fluxes
figure
plot(Qa, 'LineWidth', 2)
hold on
% plot(-Frict, 'LineWidth', 2);
% 
% plot(-Advt, 'LineWidth', 2);
% % plot(-qdivt, 'LineWidth', 2);
% 
% plot(-Diat, 'LineWidth', 2);

% plot(-(Fricst+Diast), 'LineWidth', 3, 'LineStyle', '--');
% plot(qdira);
plot(Qda, 'LineWidth', 2, 'LineStyle', '--')
plot(-Fricst); 
plot(-Diast);
plot(-(Fricst+Diast));

plot(-(Frict+Diat));
plot(-FricstN);
plot(-DiastN);
plot(-(FricstN+DiastN));
% plot(-(Qa-Qda));
hold off
legend('Q','Q_{Direct}',  'Fric(0)', 'Dia(0)', 'Sum(0)', 'Sum(\nabla J)', 'FricN', 'DiaN');
xlabel('Num time steps (60 min)');
ylabel('\Delta Q');
grid on

%%
% figure
% plot((Qa+Fricst+Diast)./Qa)
% scatter(Qa, -Fricst-Diast)
%%
% Make Time Series Figure of DeltaQ and Fluxes
cl = [15.8 16.85];
time = T(2:end-1)./86400; %in days
fs = 16;
xpos = 16;

figure
subplot(2,2,1:2)
plot(time, Qa, 'LineWidth', 3)
hold on
plot(time, -Fricst, 'LineWidth', 2); 
plot(time, -Diast, 'LineWidth', 2);
plot(time, -(Fricst+Diast), 'LineWidth',3, 'LineStyle', '--');
% plot(-(Fricst+Diast+FricstN + DiastN))
% plot(-(Frict + Diat));
% plot(Qda+Advt);
hold off
legend('Q', '-J_F^z', '-J_B^z','Sum', 'Location', 'SouthWest');
xlabel('Days');
ylabel('\Delta Q');
grid on
title(['Surface B_0: ', num2str(squeeze(B0(1,1,1))),'    Temperature Range: ', num2str(isoT(1)),'-',num2str(isoT(2))])
set(gca, 'FontSize', fs);

subplot(2,2,3)
contourf(X/1000, Y/1000, squeeze(THETA(:,:,1,tind)).'); shading interp
set(gca, 'clim', cl)
xlabel('x (km)'); ylabel('y (km)');
hold on
contour(squeeze(THETA(:,:,1,tind)).', isoT, 'k')

axis equal
grid on
set(gca, 'xlim', [0 X(end)/1000], 'ylim', [0 Y(end)/1000]);
set(gca, 'FontSize', fs);
yt = Y/1000;
plot(xpos*ones(size(yt)), yt,'--k', 'LineWidth', 2)
hold off
title(['Day: ', num2str(time(tind), 2)])

subplot(2,2,4)
[c, h]=contourf(Y/1000, Z, squeeze(THETA(xpos,:,:,tind)).', 20); 
% set(h, 'edgecolor', 'none');
set(gca, 'clim', cl)
xlabel('y (km)'); ylabel('z (m)');
hold on
contour(Y/1000, Z, squeeze(THETA(xpos,:,:,tind)).', isoT, 'k')
% set(gca, 'ydir', 'reverse')
cb = colorbar;
set(get(cb, 'ylabel'), 'string', 'T (\circ C)');
set(gca, 'FontSize', fs);
grid on

set(gcf, 'Color', 'w', 'Position', [675   370   825   605]);





%%
figure
subplot(3,1,1)
% pcolor(time, Y/1000, squeeze(THETA(xpos,:,1,:))); shading interp
[c h] = contourf(time, Y/1000, squeeze(THETA(xpos,:,1,:)), 50); 
set(h, 'edgecolor', 'none')

ylabel('y (km)');
xlabel('Days');
% colorbar;
grid on
set(gca, 'FontSize', fs);

subplot(3,1,2:3)

plot(time, smooth(Qt,1), 'LineWidth', 2)
hold on
plot(time, -Frics, 'LineWidth', 2);
% plot(-Adv, 'LineWidth', 2);
plot(time, -Dias, 'LineWidth', 2);


plot(time,-(Frics+Dias), 'LineWidth', 2, 'LineStyle', '--');
% plot(qdira);
% plot(Qta, 'LineWidth', 2)

hold off
legend('\partialQ/\partial t', 'Fric', 'Dia', 'Sum',  'Location', 'SouthEast');
xlabel('Days');
ylabel('\partialQ/\partial t');
grid on
set(gca, 'FontSize', fs);
% cb = colorbar;
% set(cb, 'visible', 'off')
set(gcf, 'Color', 'w', 'Position', [675   445   830   530]);

%%
% Pointwise budget
[~, ~, ~, qt] = gradient(Q, ts);
res = qt - FricDiv - AdvDiv - DiaDiv;
% res = qt - Qdir;
%%
indj = 24;
indk = 46;
indd = 3;

figure
plot(squeeze(qt(indj, indk, indd,:)), 'LineWidth', 2);
hold on
% indd=10;
plot(squeeze(-AdvDiv(indj, indk, indd,:) - FricDiv(indj, indk, indd,:) - DiaDiv(indj, indk, indd,:)), 'LineWidth',2, 'LineStyle', '--');
plot(squeeze(-AdvDiv(indj, indk, indd,:)) );
plot(squeeze( - FricDiv(indj, indk, indd,:) ));
plot(squeeze( - DiaDiv(indj, indk, indd,:)));

hold off

legend('dqdt', 'RHS SUM', '-Div(J_A)', '-Div(J_F)', '-Div(J_B)');

%%
divs = -AdvDiv -DiaDiv-FricDiv;
figure
subplot(1,2,1)
plot(squeeze(nanmean(nanmean(nanmean(qt,4),3))))
hold on
plot(squeeze(nanmean(nanmean(nanmean(divs,4),3))));
plot(squeeze(nanmean(nanmean(nanmean(-FricDiv,4),3))));
plot(squeeze(nanmean(nanmean(nanmean(-DiaDiv,4),3))));


hold off

subplot(1,2,2)
pcolor(squeeze(nanmean(nanmean(qt-divs,3)))); shading interp
hold on
% pcolor(squeeze(nanmean(nanmean(divs,3))));
hold off


toc