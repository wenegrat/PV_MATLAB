%%
THETA = outputFull.T;
[nx, ny,nz, nt] = size(THETA);
Z = outputFull.Z;
Zfull = permute(repmat(Z, [1,nx, ny, nt]), [2 3 1 4]);

%%
U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);
%%
KE = calculateKE(U, V, W, X, Y, Z, time);
[APE, PE, SC] = calculateAPE(THETA, X, Y, Z, time);


%%
% APEm =squeeze( nanmean(nanmean(nanmean(APE))));

%%
xl = [0 30];

gap = [0.05 0.01]; margh = 0.1; margw=0.1;
figure
subtightplot(2,1,1, gap, margh, margw)
plot(time, output.dJba_t, 'LineWidth', 2);
hold on
plot(time, output.dJbsa, 'LineWidth', 2);
plot(time, output.dJbea,'LineWidth', 2);
hold off
l=legend('$J_{D_{Total}}$', '$J_{D_{SURF}}$', '$J_{D_{EDDY}}$');
set(l, 'Interpreter', 'Latex', 'FontSize', 20);
grid on

set(gca, 'xlim', xl);
set(gca, 'XTickLabel', []);
set(gca, 'FontSize', fs)
ylabel('$m$ $s^{-4}$', 'FontSize', 19);
subtightplot(2,1,2, gap, margh, margw)
KEd = -(KE-KE(1));
APEd = -(APE.'-APE(1));
plot(time, gradient(KEd, ts),'LineWidth', 2);
hold on
plot(time, gradient(APEd, ts),'LineWidth', 2);
% plot(time, -(KE-KE(1))-(APE.'-APE(1)));
% plot(time, PE,'LineWidth', 2);
% plot(time, cumtrapz(SC)*ts);
% plot(time, APEm-APEm(1) + (KEm-KEm(1)))
hold off
grid on
l = legend('$\frac{\partial KE}{\partial t}$', '$\frac{\partial APE}{\partial t}$');
set(l, 'Interpreter', 'Latex', 'FontSize', 20);
set(gca, 'xlim', xl);
set(gca, 'FontSize', fs)
xlabel('Days', 'FontSize', 20);
ylabel('$kg$ $m^{-1}s^{-3}$', 'FontSize', 19);
set(gcf, 'Color', 'w', 'Position', [ 675   478   716   496])