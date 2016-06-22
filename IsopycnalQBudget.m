%% 
% GENERATE ISOPYCNAL BUDGET
  clc;
 statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
% Parameters            
dx = 1000; dy = dx; dz = 2.5;
nx = 48; ny=96; nz=200;
ts = 1800;
TtoB = 9.81.*2e-4;
tslice = [300 1300];
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
mask = zeros(nx, ny, nz, tslice(end)-tslice(1)+1);
isoT = [16.55 16.8];
isoT = [16.25 16.55];
isoT = [0 30];
for i=1:(tslice(end)-tslice(1)+1);
%    slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};
   mask(:,:,:,i) = (THETA(:,:,:,i)>isoT(1)) & (THETA(:,:,:,i)<isoT(2));
end 

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
withsides = false;
Fi = JFz(:,:,zl,:).*mask(:,:,zl,:);

Fi(~isfinite(Fi)) = 0;
% Fi = cumsum(Fi, 4).*ts;
Frics = squeeze(nansum(nansum(JFz(:,yl:end-yl,zl,:).*mask(:,yl:end-yl,zl,:))) - nansum(nansum(JFz(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dx.*dy;
%Frics = squeeze(nansum(nansum(JFz(:,yl:end-yl,zl,:).*mask(:,yl:end-yl,zl,:))) ).*dx.*dy;

% Frics = squeeze(trapz(trapz(Fi))).*dx.*dy;
if withsides; Frics = Frics+squeeze(nansum(nansum(JFy(:,end-yl,:,:).*mask(:,end-yl,:,:))) - nansum(nansum(JFy(:,yl,:,:).*mask(:,yl,:,:)))).*dz.*dx; end;
% Frics = Frics./vol;
% Fricst = Frics./vol;
Fricst = cumtrapz(Frics).*ts./vol;
Fricst = Fricst - Fricst(1);
Di = JBz(:,:,zl,:).*mask(:,:,zl,:);
Di(~isfinite(Di))=0;
% Dias = squeeze(trapz(trapz(Di))).*dx.*dy;
Dias = squeeze(nansum(nansum(JBz(:,yl:end-yl,zl,:).*mask(:,yl:end-yl,zl,:))) - nansum(nansum(JBz(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dx.*dy;
%Dias = squeeze(nansum(nansum(JBz(:,yl:end-yl,zl,:).*mask(:,yl:end-yl,zl,:)))).*dx.*dy;

if withsides; Dias = Dias+squeeze(nansum(nansum(JBy(:,end-yl,:,:).*mask(:,end-yl,:,:))) - nansum(nansum(JBy(:,yl,:,:).*mask(:,yl,:,:)))).*dz.*dx; end;
% Dias = Dias./vol;
Diast = cumtrapz(Dias).*ts./vol;
Diast = Diast - Diast(1);

Qi = Qdir.*mask;
Qta = squeeze(nansum(nansum(nansum(Qi)))).*gridvol; %This is volume integral of dQ/dt
Qda = cumtrapz(Qta).*ts./vol; %Time Integral of Vol Integral of dQ/dt.
Qda = Qda - Qda(1);

Qi = Q.*mask;
% Qi(:,:,1,:) = Qi(:,:,2,:);
% Qi(~isfinite(Qi)) = 0;
% Qi(Qi>5e-9) = 0; % XXX-HACK TEST
% Qt = gradient`
% Qa = squeeze(trapz(trapz(trapz(Qi)))).*gridvol; %Note that doing it this way is incorrect, as it doesn't account for time variations in Vol...
Qa = squeeze(nansum(nansum(nansum(Qi)))).*gridvol;
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
plot(-Fricst); 
plot(-Diast);
plot(-(Fricst+Diast));
% plot(Qda+Advt);
hold off
legend('Q', 'Fric', 'Dia', 'Sum');
xlabel('Num time steps (Hourly)');
ylabel('\Delta Q');
grid on

%%
% Make Time Series Figure of DeltaQ and Fluxes
figure
plot(Qa, 'LineWidth', 2)
hold on
plot(-Frict, 'LineWidth', 2);

plot(-Advt, 'LineWidth', 2);
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
plot(smooth(Qt, 1), 'LineWidth', 2)
hold on
plot(-Frics, 'LineWidth', 2);
% plot(-Adv, 'LineWidth', 2);
plot(-Dias, 'LineWidth', 2);


plot(-(Frics+Dias), 'LineWidth', 2, 'LineStyle', '--');
% plot(qdira);
plot(Qta, 'LineWidth', 2)

hold off
legend('Q', 'Fric', 'Dia', 'Sum', 'Qdir');
xlabel('Num time steps (Hourly)');
ylabel('\Delta Q');
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