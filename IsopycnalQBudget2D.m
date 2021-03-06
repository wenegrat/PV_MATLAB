%% 
% GENERATE ISOPYCNAL BUDGET
  clc;
 statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
% Parameters            
dx = 2500; dy = dx; dz = 2.5;
nx = 1; ny=20; nz=200;
ts = 3600;
TtoB = 9.81.*2e-4;
tslice = [10 200];
tslice = [10 479];
slice={0, 0, 0, tslice};%100 120
sliceEta={0,0,[1 1],tslice};%251 271

ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, nx, ny, 1]), [2 3 1 4]);
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

for i=1:(tslice(end)-tslice(1)+1)
    disp(num2str(i))
   slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};
   [Q(:,:,:,i), Qdir(:,:,:,i), JAx(:,:,:,i), JAy(:,:,:,i), JAz(:,:,:,i), ...
       JFx(:,:,:,i), JFy(:,:,:,i), JFz(:,:,:,i), JBx(:,:,:,i), JBy(:,:,:,i), JBz(:,:,:,i)] ...
       = calcQBudget2D(diagfile, statefile, etanfile, [nx, ny], slicetemp, dx, dy);

   
   FricDiv(:,:,:,i) = + Drv(dy, JFy(:,:,:,i), 'y') + Drv(metric, JFz(:,:,:,i), 'z');
   AdvDiv(:,:,:,i) = + Drv(dy, JAy(:,:,:,i), 'y') + Drv(metric, JAz(:,:,:,i), 'z');
   DiaDiv(:,:,:,i) =  + Drv(dy, JBy(:,:,:,i), 'y') + Drv(metric, JBz(:,:,:,i), 'z');
   THETA(:,:,:,i) = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);

end

%%
% Generate Mask
mask = zeros(nx, ny, nz, tslice(end)-tslice(1)+1);
isoT = [16.55 16.9];
isoT = [16.5 16.95];
 isoT = [3 30];
for i=1:(tslice(end)-tslice(1)+1);
%    slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};
   mask(:,:,:,i) = (THETA(:,:,:,i)>isoT(1)) & (THETA(:,:,:,i)<isoT(2));
end 
% 
% [c, h] = contourf(squeeze(THETA(:,:,1,700)).', linspace(15, 20, 100)); shading interp
% set(h, 'edgecolor','none')
% xlabel('x'); ylabel('y');
% hold on
% contour(squeeze(THETA(:,:,1,700)).', isoT, 'k')
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



% JAxi = JAx.*mask; JAyi = JAy.*mask; JAzi = JAz.*mask;
% JFxi = JFx.*mask; JFyi = JFy.*mask; JFzi = JFz.*mask;
% JBxi = JBx.*mask; JByi = JBy.*mask; JBzi = JBz.*mask;
Frici = FricDiv.*mask;
Advi = AdvDiv.*mask;
Diai = DiaDiv.*mask;

vol = dx.*dy.*dz.*squeeze(sum(sum(sum(mask))));
% vol = 1;
gridvol = dx.*dy.*dz;
Fric = squeeze(nansum(nansum(nansum(Frici)))).*gridvol;
Frict = cumtrapz(Fric).*ts./vol;
Adv = squeeze(nansum(nansum(nansum(Advi)))).*gridvol;
Advt = cumtrapz(Adv).*ts./vol;
Dia = squeeze(nansum(nansum(nansum(Diai)))).*gridvol;
Diat = cumtrapz(Dia).*ts./vol;

zl =1;
yl = 2;
withsides = true;
%Original Flux Calc
Frics = squeeze((nansum(JFz(:,yl:end-yl,zl,:).*mask(:,yl:end-yl,zl,:))) - (nansum(JFz(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dx.*dy;
if withsides; Frics = Frics+squeeze((nansum(JFy(:,end-yl,:,:).*mask(:,end-yl,:,:))) -(nansum(JFy(:,yl,:,:).*mask(:,yl,:,:)))).*dz.*dx; end;
Fricst = cumtrapz(Frics).*ts./vol;
Fricst = Fricst - Fricst(1);

Dias = squeeze((nansum(JBz(:,yl:end-yl,zl,:).*mask(:,yl:end-yl,zl,:))) - (nansum(JBz(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dx.*dy;
if withsides; Dias = Dias+squeeze((nansum(JBy(:,end-yl,:,:).*mask(:,end-yl,:,:))) - (nansum(JBy(:,yl,:,:).*mask(:,yl,:,:)))).*dz.*dx; end;
Diast = cumtrapz(Dias).*ts./vol;
Diast = Diast - Diast(1);

%Alternate Flux Calc
% t = (tslice(1):1:tslice(end))*ts;
% JINT = cumtrapz(t, JFz, 4);
% Frics = squeeze(nansum(nansum(JINT(:,yl:end-yl, zl,:).*mask(:,yl:end-yl,zl,:))) - nansum(nansum(JINT(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dy.*dy;
% Fricst = Frics./vol;
% 
% BINT = cumtrapz(t, JBz, 4);
% Dias = squeeze(nansum(nansum(BINT(:,yl:end-yl, zl,:).*mask(:,yl:end-yl,zl,:))) - nansum(nansum(BINT(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dy.*dy;
% Diast = Dias./vol;
% 
Qi = Qdir.*mask;
Qta = squeeze(nansum(nansum(nansum(Qi)))).*gridvol; %This is volume integral of dQ/dt
Qda = cumtrapz(Qta).*ts./vol; %Time Integral of Vol Integral of dQ/dt.
Qda = Qda - Qda(1);

Qi = Q.*mask;
% Qi(:,:,1,:) = Qi(:,:,2,:);
Qi(~isfinite(Qi)) = 0;
% Qi(Qi>5e-9) = 0; % XXX-HACK TEST
% Qt = gradient`
% Qa = squeeze(trapz(trapz(trapz(Qi)))).*gridvol; %Note that doing it this way is incorrect, as it doesn't account for time variations in Vol...
Qa = squeeze((nansum(nansum(Qi)))).*gridvol;
% Qda = cumtrapz(Qta).*ts; %Need to check if this is correct.
Qt = gradient(Qa, ts);
Qa = (Qa - Qa(1))./vol;% 
% Jbah = squeeze(nansum(nansum(JBz(:,ylims,:,:))))*dx*dy./area;
% Jfah = squeeze(nansum(nansum(JFz(:,ylims,:,:))))*dx*dy./area;
% qah = squeeze(nansum(nansum(Q(:,ylims,:,:)))).*dx.*dy./area;
% % qah(:,:,1:2,:) = NaN;
% % qdirah = squeeze(nansum(nansum(Qdir))).*dx.*dy./area;
% 
% 
% %Meridional Components
% Jfahy = squeeze(nansum(JFy(:, ylims(end), :,:)-JFy(:,ylims(1),:,:))).*dx./area;
% Jfay = nansum(Jfahy(zlims,:)).*dz;
% Jfay = cumtrapz(Jfay).*ts;
% Jaahy = squeeze(nansum(JAy(:, ylims(end), :,:)-JAy(:,ylims(1),:,:))).*dx./area;
% Jaay = nansum(Jaahy(zlims,:)).*dz;
% Jaay = cumtrapz(Jaay).*ts;
% 
% %Integrating the vertical difference in time.
% deltaJb = Jbah(zlims(1),:) - Jbah(zlims(end),:);
% % deltaJb = Jbah(2,:) - Jbah(end,:);
% 
% Jbza = cumtrapz(deltaJb).*ts;
% deltaJf = Jfah(1,:) - Jfah(end,:);
% deltaJf = Jfah(zlims(1),:) - Jfah(zlims(end),:);
% 
% Jfza = cumtrapz(deltaJf).*ts;

% Jfza = cumtrapz(Jfah(2,:)).*ts;

%Integrate vertically`
% qa = nansum(qah(zlims,:))*dz;
% qt = gradient(qa, ts);
% 
% qa = qa - qa(1); %Define delta quantity


% qa = -nansum(qah)*dz; %integrate in z
% qa = cumtrapz(qa)*ts; %integrate in time.

% qdira = -nansum(qdirah)*dz; %integrate in 
% qdira = cumtrapz(qdira)*ts; %integrate in time.
%%
% Make Time Series Figure of DeltaQ and Fluxes
figure
subplot(1,4,1:3)
plot(Qa, 'LineWidth', 2)
hold on
plot(-Fricst); 
plot(-Diast);
plot(-(Fricst+Diast));
% plot(Qda+Advt);
hold off
legend('Q', 'Fric', 'Dia', 'Sum');
xlabel('Num time steps (Hours)');
ylabel('\Delta Q');
grid on
title(num2str(isoT))
% pcolor(squeeze(THETA(:,:,1,700)).'); shading interp
% set(gca, 'clim', [16 17])
% xlabel('x'); ylabel('y');
% hold on
% contour(squeeze(THETA(:,:,1,700)).', isoT, 'k')
% colorbar
[nx, ny, nz, nt] = size(Q);
subplot(1,4,4)
pcolor(squeeze(THETA(1,:,:,floor(nt))).'); shading interp
set(gca, 'clim', [16 17])
xlabel('y'); ylabel('z');
hold on
contour(squeeze(THETA(1,:,:,floor(nt))).', isoT, 'k')
set(gca, 'ydir', 'reverse')
colorbar
set(gcf, 'Position', [355          82        1479         533]);

%%
% Normalized Error Plot
figure
err = Qa + (Fricst+Diast);

plot(err./Qa, 'LineWidth',2)
title('Res/Q');
xlabel('Time Step');
grid on
%%
% Make Time Series Figure of DeltaQ and Fluxes
figure
plot(Qa, 'LineWidth', 2)
hold on
plot(-Frict, 'LineWidth', 2);

% plot(-Advt, 'LineWidth', 2);
% plot(-qdivt, 'LineWidth', 2);

plot(-Diat, 'LineWidth', 2);

plot(-(Frict+Diat), 'LineWidth', 3, 'LineStyle', '--');
% plot(qdira);
plot(Qda, 'LineWidth', 2)
plot(-Fricst); 
plot(-Diast);
plot(-(Fricst+Diast));
plot(-(Qa-Qda));
hold off
legend('Q', 'Fric','Adv', 'Dia', 'Sum (w adv)', 'Qdir', 'Fric(0)', 'Dia(0)', 'Sum(0)');
xlabel('Num time steps (Hourly)');
ylabel('\Delta Q');
grid on



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

scatter(Qa, -Fricst-Diast)
% hold on
% scatter(Qa+(Fricst+Diast), (vdivat).*gridvol./vol)
% hold off
grid on

%%
plot(Qta);
hold on
plot(Qt./vol); hold off
grid on