%% 
% GENERATE ISOPYCNAL BUDGET
  clc; clear all; close all;
 statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
% Parameters            
% dx = 500; dy = dx; dz = 3;
% nx = 96; ny=192; nz=200;
% ts = 3600;

dx = 1000; dy = dx; dz = 2.5;
nx = 48; ny=96; nz=200;
ts = 1800;

TtoB = 9.81.*2e-4;
tslice = [10 200];
tslice = [10 69];
tslice = [600 699];
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
mask = Q;
FricDiv = Q;
AdvDiv = Q;
DiaDiv = Q;
THETA  = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);
resid = Q;
uterm = Q;
bterm = Q;
% remaining = (tslice(end)-tslice(1)+1);
incc =10;
inc = incc-1;

ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, nx, ny, incc]), [2 3 1 4]);
for i=1:incc:(tslice(end)-tslice(1)+1)
    disp(num2str(i));
%     remaining = @(x)    remaining -1
%     disp(['Remaining: ', num2str(remaining), '    (Processing: ', num2str(i),')'])
   slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1+inc]};
   [Q(:,:,:,i:i+inc), Qdir(:,:,:,i:i+inc), JAx(:,:,:,i:i+inc), JAy(:,:,:,i:i+inc), JAz(:,:,:,i:i+inc), ...
       JFx(:,:,:,i:i+inc), JFy(:,:,:,i:i+inc), JFz(:,:,:,i:i+inc), JBx(:,:,:,i:i+inc), JBy(:,:,:,i:i+inc), JBz(:,:,:,i:i+inc)] ...
       = calcQBudget(diagfile, statefile, etanfile, [nx, ny, incc], slicetemp, dx, dy);

   FricDiv(:,:,:,i:i+inc) = Drv(dx, JFx(:,:,:,i:i+inc), 'x') + Drv(dy, JFy(:,:,:,i:i+inc), 'y') + Drv(metric, JFz(:,:,:,i:i+inc), 'z');
   AdvDiv(:,:,:,i:i+inc) = Drv(dx, JAx(:,:,:,i:i+inc), 'x') + Drv(dy, JAy(:,:,:,i:i+inc), 'y') + Drv(metric, JAz(:,:,:,i:i+inc), 'z');
   DiaDiv(:,:,:,i:i+inc) = Drv(dx, JBx(:,:,:,i:i+inc), 'x') + Drv(dy, JBy(:,:,:,i:i+inc), 'y') + Drv(metric, JBz(:,:,:,i:i+inc), 'z');
   THETA(:,:,:,i:i+inc) = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);
   [resid(:,:,:,i:i+inc), uterm(:,:,:,i:i+inc), bterm(:,:,:,i:i+inc)] = AssessQCancellation(statefile, diagfile, slicetemp);
end

%%
% Generate Mask
mask = zeros(nx, ny, nz, tslice(end)-tslice(1)+1);
isoT = [16.55 16.9];
% isoT = [16.5 16.95];
 isoT = [3 30];
% isoT = [16.9 30];
% isoT = [3 16.55];
[nx, ny, nz, nt] = size(THETA);
for i=1:nt;
%    slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};
   mask(:,:,:,i) = (THETA(:,:,:,i)>isoT(1)) & (THETA(:,:,:,i)<isoT(2));
end 

mask(:,end-2:end,:,:) = 0; %wall effects that might be in Q
%%
tind = floor((tslice(2)-tslice(1))/2);

[c, h] = contourf(squeeze(THETA(:,:,1,tind)).', linspace(15, 20, 100)); shading interp
set(h, 'edgecolor','none')
xlabel('x'); ylabel('y');
hold on
contour(squeeze(THETA(:,:,1,tind)).', isoT, 'k')
set(gca, 'clim', [16 17])
colorbar
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
% Frici = cumtrapz(t, Frici, 4);
Frici = Frici.*mask;
% Fric = squeeze(nansum(nansum(nansum(Frici)))).*gridvol;
Fric = squeeze(trapz(trapz(trapz(Frici)))).*gridvol;
Fric = cumtrapz(t, Fric);
Frict = Fric./vol;
% Advi = cumtrapz(t, Advi, 4);
Advi = Advi.*mask;
Adv = squeeze(nansum(nansum(nansum(Advi)))).*gridvol;
Adv = cumtrapz(t,Adv);
Advt = Adv./vol;
% Diai = cumtrapz(t, Diai, 4);
Diai = Diai.*mask;
Dia = squeeze(nansum(nansum(nansum(Diai)))).*gridvol;
Fric = squeeze(trapz(trapz(trapz(Frici)))).*gridvol;

Dia = cumtrapz(t,Dia);
Diat = Dia./vol;

%Calculate as surface fluxes
zl =1;
yl = 1;
withsides = false;
%Original Flux Calc
Frics = squeeze(nansum(nansum(JFz(:,yl:end-yl,zl,:).*mask(:,yl:end-yl,zl,:))) - nansum(nansum(JFz(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dx.*dy;
if withsides; Frics = Frics+squeeze(nansum(nansum(JFy(:,end-yl,:,:).*mask(:,end-yl,:,:))) - nansum(nansum(JFy(:,yl,:,:).*mask(:,yl,:,:)))).*dz.*dx; end;
Fricst = cumtrapz(Frics).*ts./vol;
Fricst = Fricst - Fricst(1);

Dias = squeeze(nansum(nansum(JBz(:,yl:end-yl,zl,:).*mask(:,yl:end-yl,zl,:))) - nansum(nansum(JBz(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dx.*dy;
if withsides; Dias = Dias+squeeze(nansum(nansum(JBy(:,end-yl,:,:).*mask(:,end-yl,:,:))) - nansum(nansum(JBy(:,yl,:,:).*mask(:,yl,:,:)))).*dz.*dx; end;
Diast = cumtrapz(Dias).*ts./vol;
Diast = Diast - Diast(1);

%Alternate Flux Calc
% JINT = cumtrapz(t, JFz, 4);
% % Frics = squeeze(nansum(nansum(JINT(:,yl:end-yl, zl,:).*mask(:,yl:end-yl,zl,:))) - nansum(nansum(JINT(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dy.*dy;
% Frics = squeeze(nansum(nansum(JINT(:,yl:end-yl, zl,:).*mask(:,yl:end-yl,zl,:))) ).*dy.*dy;
% 
% Fricst = Frics./vol;
% 
% BINT = cumtrapz(t, JBz, 4);
% % Dias = squeeze(nansum(nansum(BINT(:,yl:end-yl, zl,:).*mask(:,yl:end-yl,zl,:))) - nansum(nansum(BINT(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dy.*dy;
% Dias = squeeze(nansum(nansum(BINT(:,yl:end-yl, zl,:).*mask(:,yl:end-yl,zl,:))) ).*dy.*dy;
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
Qi(~isfinite(Qi)) = 0;
Qi = Qi.*mask;
Qa = squeeze(nansum(nansum(nansum(Qi)))).*gridvol;
Qt = gradient(Qa, ts);
Qa = Qa-Qa(1);
Qa = Qa./vol;% 
% Qa = Qa./1035
%%
residual = cumtrapz(t, uterm+bterm, 4);
residual = squeeze(nansum(nansum(nansum(residual.*mask)))).*gridvol;
residual = residual./vol;
residavg = residual;
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

plot(-(Fricst+Diast), 'LineWidth', 3, 'LineStyle', '--');
% plot(qdira);
plot(Qda, 'LineWidth', 2, 'LineStyle', '--')
plot(-Fricst); 
plot(-Diast);
plot(-(Fricst+Diast));

plot(-(Frict+Diat));
% plot(-(Qa-Qda));
hold off
legend('Q', 'Sum (no Res)','Q_{Direct}',  'Fric(0)', 'Dia(0)', 'Sum(0)', 'Sum(\nabla J)');
xlabel('Num time steps (30 min)');
ylabel('\Delta Q');
grid on

%%
% scatter(Qa, -Fricst-Diast)
%%
% Make Time Series Figure of DeltaQ and Fluxes
figure
subplot(2,2,1:2)
plot(Qa, 'LineWidth', 2)
hold on
plot(-Fricst); 
plot(-Diast);
plot(-(Fricst+Diast));
% plot(Qda+Advt);
hold off
legend('Q', 'Fric', 'Dia', 'Sum');
xlabel('Num time steps (30 minutes)');
ylabel('\Delta Q');
grid on
title(num2str(isoT))
subplot(2,2,3)
pcolor(squeeze(THETA(:,:,1,tind)).'); shading interp
set(gca, 'clim', [16 17])
xlabel('x'); ylabel('y');
hold on
contour(squeeze(THETA(:,:,1,tind)).', isoT, 'k')
colorbar
subplot(2,2,4)
pcolor(squeeze(THETA(20,:,:,tind)).'); shading interp
set(gca, 'clim', [16 17])
xlabel('x'); ylabel('z');
hold on
contour(squeeze(THETA(20,:,:,tind)).', isoT, 'k')
set(gca, 'ydir', 'reverse')
colorbar





%%
figure
plot(smooth(Qt,1), 'LineWidth', 2)
hold on
plot(-Frics, 'LineWidth', 2);
% plot(-Adv, 'LineWidth', 2);
plot(-Dias, 'LineWidth', 2);


plot(-(Frics+Dias), 'LineWidth', 2, 'LineStyle', '--');
% plot(qdira);
% plot(Qta, 'LineWidth', 2)

hold off
legend('dQ/dt', 'Fric', 'Dia', 'Sum', 'Qdir');
xlabel('Num time steps (30 min)');
ylabel('dQ/dt');
grid on

%%
% Pointwise budget
[~, ~, ~, qt] = gradient(Q, ts);
%%
indj = 45;
indk = 45;
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
