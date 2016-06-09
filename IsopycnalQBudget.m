%% 
% GENERATE ISOPYCNAL BUDGET
  clc;
 statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
% Parameters            
dx = 1000; dy = dx; dz = 2.5;
nx = 96; ny=nx; nz=200;
ts = 1800;
TtoB = 9.81.*2e-4;
tslice = [1 1339];
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
       = calcQBudget(diagfile, statefile, etanfile, [nx, ny], slicetemp, dx, dy);

   
   FricDiv(:,:,:,i) = Drv(dx, JFx(:,:,:,i), 'x') + Drv(dy, JFy(:,:,:,i), 'y') + Drv(metric, JFz(:,:,:,i), 'z');
   AdvDiv(:,:,:,i) = Drv(dx, JAx(:,:,:,i), 'x') + Drv(dy, JAy(:,:,:,i), 'y') + Drv(metric, JAz(:,:,:,i), 'z');
   DiaDiv(:,:,:,i) = Drv(dx, JBx(:,:,:,i), 'x') + Drv(dy, JBy(:,:,:,i), 'y') + Drv(metric, JBz(:,:,:,i), 'z');
   THETA(:,:,:,i) = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);

end

%%
% Generate Mask
isoT = [16.5 16.55];
% isoT = [0 20];
for i=1:(tslice(end)-tslice(1)+1);
%    slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};
   mask(:,:,:,i) = (THETA(:,:,:,i)>isoT(1)) & (THETA(:,:,:,i)<isoT(2));
end
%%
% Calculate divergences
% Integrate in Volume.
disp('Integrate Over Volume');

%Average Horizontally
ylims = 2:(ny-2);
zlims = 2:nz-1;


% JAxi = JAx.*mask; JAyi = JAy.*mask; JAzi = JAz.*mask;
% JFxi = JFx.*mask; JFyi = JFy.*mask; JFzi = JFz.*mask;
% JBxi = JBx.*mask; JByi = JBy.*mask; JBzi = JBz.*mask;
Frici = FricDiv.*mask;
Advi = AdvDiv.*mask;
Diai = DiaDiv.*mask;

vol = dx.*dy.*dz.*squeeze(sum(sum(sum(mask))));
gridvol = dx.*dy.*dz;
Fric = squeeze(nansum(nansum(nansum(Frici)))).*gridvol./vol;
Frict = cumtrapz(Fric).*ts;
Adv = squeeze(nansum(nansum(nansum(Advi)))).*gridvol./vol;
Advt = cumtrapz(Adv).*ts;
Dia = squeeze(nansum(nansum(nansum(Diai)))).*gridvol./vol;
Diat = cumtrapz(Dia).*ts;

zl  = 1;
yl = 2;
Frics = squeeze(nansum(nansum(JFz(:,:,zl,:).*mask(:,:,zl,:))) - nansum(nansum(JFz(:,:,end-zl+1,:).*mask(:,:,end-zl+1,:)))).*dx.*dy;
Frics = Frics+squeeze(nansum(nansum(JFz(:,end-yl,:,:).*mask(:,end-yl,:,:))) - nansum(nansum(JFz(:,yl,:,:).*mask(:,yl,:,:)))).*dz.*dx;
Frics = Frics./vol;
Fricst = cumtrapz(Frics).*ts;
% Dias = squeeze(nansum(nansum(JBz(:,:,zl,:).*mask(:,:,zl,:))));
Dias = squeeze(nansum(nansum(JBz(:,:,zl,:).*mask(:,:,zl,:))) - nansum(nansum(JBz(:,:,end-zl+1,:).*mask(:,:,end-zl+1,:)))).*dx.*dy;
Dias = Dias+squeeze(nansum(nansum(JBz(:,end-yl,:,:).*mask(:,end-yl,:,:))) - nansum(nansum(JBz(:,yl,:,:).*mask(:,yl,:,:)))).*dz.*dx;
Dias = Dias./vol;
Diast = cumtrapz(Dias).*ts;

Qi = Qdir.*mask;
Qta = squeeze(nansum(nansum(nansum(Qi)))).*gridvol./vol; %This is volume integral of dQ/dt
Qda = cumtrapz(Qta).*ts; %Need to check if this is correct.
Qda = Qda - Qda(1);

Qi = Q.*mask;
% Qt = gradient
Qa = squeeze(nansum(nansum(nansum(Qi)))).*gridvol./vol; %Note that doing it this way is incorrect, as it doesn't account for time variations in Vol...

% Qda = cumtrapz(Qta).*ts; %Need to check if this is correct.
Qt = gradient(Qa, ts);
Qa = Qa - Qa(1);% 
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

%Integrate vertically
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
plot(Qa, 'LineWidth', 2)
hold on
plot(-Frict, 'LineWidth', 2);

plot(-Advt, 'LineWidth', 2);
plot(-Diat, 'LineWidth', 2);

plot(-(Frict+Diat), 'LineWidth', 3, 'LineStyle', '--');
% plot(qdira);
plot(Qda, 'LineWidth', 2)
plot(-Fricst); 
plot(-Diast);
plot(-(Fricst+Diast));
hold off
legend('Q', 'Fric','Adv', 'Dia', 'Sum', 'Qdir', 'Fric(0)', 'Dia(0)', 'Sum(0)');
xlabel('Num time steps (Hourly)');
ylabel('\Delta Q');
grid on

%%
figure
plot(Qt, 'LineWidth', 2)
hold on
plot(-Fric, 'LineWidth', 2);
plot(-Adv, 'LineWidth', 2);
plot(-Dia, 'LineWidth', 2);

plot(-(Fric+Adv+Dia), 'LineWidth', 2, 'LineStyle', '--');
% plot(qdira);
plot(Qta, 'LineWidth', 2)

hold off
legend('Q', 'Fric','Adv', 'Dia', 'Sum');
xlabel('Num time steps (Hourly)');
ylabel('\Delta Q');
grid on


% %%
% plot(Qa);
% plot(