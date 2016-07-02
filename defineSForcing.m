ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);

U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);
b = GetVar(statefile, diagfile, {'b', '(1)'}, slice);

OMEGAX = GetVar(statefile, diagfile, {'VVEL', ' - Dz(1)'}, slice);
OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);
OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL', '(1) + Dx(2) - Dy(3)'}, slice);

VADVTERM = GetVar(statefile, diagfile, { 'Vm_Advec', '-(1)'}, slice);
UADVTERM = GetVar(statefile, diagfile, {'Um_Advec', '-(1)'}, slice);
QXz = Drv(metric, UADVTERM, 'z',1,1);
QXy = Drv(dy, UADVTERM, 'y', 1,1);
QVz = Drv(metric, VADVTERM, 'z', 1, 1);
QVx = Drv(dx, VADVTERM, 'x',1,1);
Qmom = -bx.*QVz + by.*QXz + bz.*(QVx - QXy);

LHSb = TtoB.*GetVar(statefile, diagfile, { 'UDIAG1', '-(1)'}, slice);
LHSb = U.*bx + V.*by + W.*bz;
LHBx = Drv(dx, LHSb, 'x', 1, 1);
LHBy = Drv(dy, LHSb, 'y', 1, 1);
LHBz = Drv(metric, LHSb, 'z', 1, 1);
Qb = OMEGAX.*LHBx + OMEGAY.*LHBy + OMEGAZ.*LHBz;
Qpartial = Qmom + Qb;


Qtemp = OMEGAX.*bx + OMEGAY.*by + OMEGAZ.*bz;
Qdiver = U.*Drv(dx, Qtemp, 'x') + V.*Drv(dy, Qtemp, 'y') + W.*Drv(metric, Qtemp, 'z');
%%
[residual, uterm, bterm] = AssessQCancellation(statefile, diagfile, slice);
%%
Sa_mom = averageQForcings(-Qmom, mask, gridvol, ts);
Sa_buoy= averageQForcings(-Qb, mask, gridvol, ts);

figure
plot(Sa_mom); 
hold on
plot(Sa_buoy)
plot(Sa_mom+Sa_buoy);
hold off
legend('\nabla b \cdot (\nabla \times (u \cdot \nabla) u)', '\omega_A \cdot \nabla(u \cdot \nabla b)', 'Sum')
grid on
title('Volume and Time integrated')
xlabel('Time steps (hours)');
%%
Sa = averageQForcings(-Qpartial, mask, gridvol, ts);
Qat = averageQForcings(Qdiver, mask, gridvol, ts);
Aat = averageQForcings(AdvDiv, mask, gridvol, ts);

res = Qa + Fricst + Diast;
RESm = -(Qpartial-AdvDiv);
RESma = Sa-Qat;
%%
scatter(Qat, Aat)
grid on
%%
figure
subplot(2,2,1)
plot(res, Sa);
hold on
% plot(res, Qat);
% plot(res, -Aat);
hold off
xlabel('$\int Q_{Model} dV + \int_0^t (\bar{J^z_F(0)} + \bar{J^z_B(0)}) dt$', 'Interpreter', 'LaTex');
ylabel('-\int\int(\nabla b \cdot (\nabla \times ADV_{Mom}) + \omega_A\cdot \nabla ADV_B)dtdV');
grid on

subplot(2,2,2)
plot(RESma);
hold on
% plot( Qat);
% plot(Aat);
hold off
grid on
xlabel('Time Step');
title('\int \int -RESm dt dV');

contlim = linspace(16 , 17, 50);
subplot(2,2,3)
pcolor(t, ztmp, squeeze(nanmean(nanmean(Qpartial(:,:,:,:)./qt))))
hold on
contour(t, ztmp, squeeze(nanmean(nanmean(THETA(:,:,:,:)))), contlim, 'k')
hold off
shading interp
set(gca, 'clim', [-2 2].*3e-16);
set(gca, 'clim', [-2 2]);

xlabel('Time Step'); ylabel('z'); 
title('Horiz Averaged dQadv/dt');

subplot(2,2,4)
pcolor(t, 1:ny, squeeze(nanmean(nanmean(Qpartial(:,:,:,:)./qt,3))))
shading interp
hold on
contour(t, 1:ny, squeeze(nanmean(nanmean(THETA(:,:,:,:), 3))), contlim, 'k')
hold off
set(gca, 'clim', [-2 2].*3e-16);
set(gca, 'clim', [-2 2]);

xlabel('Time Step'); ylabel('y');
title('Depth Averaged dQadv/dt');