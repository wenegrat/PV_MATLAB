%%
% GENERATE PV PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  clc;
 statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
% Parameters            
dx = 1000; dy = dx; dz = 5;
nx = 48; ny=nx; nz=100;
ts = 14400;
TtoB = 9.81.*2e-4;
tslice = [1 119];
slice={0, 0, 0, tslice};%100 120
sliceEta={0,0,[1 1],tslice};%251 271

%%
Q = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);
JBz = Q;
JFz = Q;
JFy = Q;
JAy = Q;
for i=1:(tslice(end)-tslice(1)+1)
    disp(num2str(i))
   slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};
   [Q(:,:,:,i), JFz(:,:,:,i), JBz(:,:,:,i), JFy(:,:,:,i), JAy(:,:,:,i)] = calcQBudget(diagfile, statefile, etanfile, [nx, ny], slicetemp);
end


%%
% Integrate in Volume.
disp('Try Volume Integral');

area = nx.*dx.*ny.*dy;
%Average Horizontally
ylims = 3:(ny-2);
zlims = 2:nz-1;
if ndims(squeeze(JFz))==4
Jbah = squeeze(nansum(nansum(JBz(:,ylims,:,:))))*dx*dy./area;
Jfah = squeeze(nansum(nansum(JFz(:,ylims,:,:))))*dx*dy./area;
qah = squeeze(nansum(nansum(Q(:,ylims,:,:)))).*dx.*dy./area;
% qah(:,:,1:2,:) = NaN;
% qdirah = squeeze(nansum(nansum(Qdir))).*dx.*dy./area;

else %2 Dimensional
    disp('Blurgh');
%     Jbah = squeeze(nansum(JBz)).*dy./(ny.*dy);
%     Jfah = squeeze(nansum(JFz)).*dy./(ny.*dy);
%     qah = squeeze(nansum(Q)).*dy./(ny.*dy);
%     qdirah  = squeeze(nansum(Qdir)).*dy./(ny*dy);
end
%Meridional Components
Jfahy = squeeze(nansum(JFy(:, ylims(end), :,:)-JFy(:,ylims(1),:,:))).*dx./area;
Jfay = nansum(Jfahy(zlims,:)).*dz;
Jfay = cumtrapz(Jfay).*ts;
Jaahy = squeeze(nansum(JAy(:, ylims(end), :,:)-JAy(:,ylims(1),:,:))).*dx./area;
Jaay = nansum(Jaahy(zlims,:)).*dz;
Jaay = cumtrapz(Jaay).*ts;

%Integrating the vertical difference in time.
deltaJb = Jbah(zlims(1),:) - Jbah(zlims(end),:);
deltaJb = Jbah(2,:) - Jbah(end,:);

Jbza = cumtrapz(deltaJb).*ts;
deltaJf = Jfah(1,:) - Jfah(end,:);
% deltaJf = Jfah(zlims(1),:) - Jfah(zlims(end),:);

Jfza = cumtrapz(deltaJf).*ts;

% Jfza = cumtrapz(Jfah(2,:)).*ts;

%Integrate vertically
qa = nansum(qah(zlims,:))*dz;
qt = gradient(qa, ts);

qa = qa - qa(1); %Define delta quantity


% qa = -nansum(qah)*dz; %integrate in z
% qa = cumtrapz(qa)*ts; %integrate in time.

% qdira = -nansum(qdirah)*dz; %integrate in z
% qdira = cumtrapz(qdira)*ts; %integrate in time.

% Make Time Series Figure of DeltaQ and Fluxes
figure
plot(qa, 'LineWidth', 2)
hold on
plot(-Jbza, 'LineWidth', 2);
plot(-Jfza, 'LineWidth', 2);

plot(-(Jbza+Jfza), 'LineWidth', 2, 'LineStyle', '--');
% plot(qdira);
plot(-Jfay);
plot(-Jaay);

hold off
legend('Q', 'J_B^z', 'J_F^z', '-SUM');
xlabel('Num time steps (Hourly)');
ylabel('\Delta Q');
grid on

%%
% In terms of dQ/dt
figure
plot(qt, 'LineWidth', 2);
hold on
plot(-deltaJb, 'LineWidth', 2);
plot(-deltaJf, 'LineWidth', 2);
plot(-(deltaJb+deltaJf), 'LineWidth',2 , 'LineStyle', '--');
hold off
grid on
legend('dQ/dt',  'J_B^z', 'J_F^z', '-SUM');