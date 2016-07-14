%% ConfirmNumHoriz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the magnitude of the numeric terms and horizontal terms relative to
% total.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isoT = [0 100]; %Needed for plot below
mask = ones(nx, ny, nz, tslice(end)-tslice(1)+1);
vol = gridvol.*squeeze(sum(sum(sum(mask))));

%INTEGRATE Q
% IntegrateQTerms;

%AreaIntegrateJTerms
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JFb, ~    ] = areaIntegrateJVecs(squeeze(JFz(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JFa = JFs-JFb;

[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JBb, ~    ] = areaIntegrateJVecs(squeeze(JBz(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JBa = JBs-JBb;

%AreaIntegrate Horizontal Terms
[JFsH, dJFHdt] = areaIntegrateJVecs(squeeze(JFzH(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JFbH, ~    ] = areaIntegrateJVecs(squeeze(JFzH(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JFaH = JFsH-JFbH;

[JBsN, ~] = areaIntegrateJVecs(squeeze(JBzN(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JBbN, ~    ] = areaIntegrateJVecs(squeeze(JBzN(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JBaN = JBsN-JBbN;

[JFsN, ~] = areaIntegrateJVecs(squeeze(JFzN(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JFbN, ~    ] = areaIntegrateJVecs(squeeze(JFzN(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JFaN = JFsN-JFbN;

%%
figure
subplot(2,1,1)
plot(time, -JFa, 'LineWidth', 2);
hold on
plot(time, -JFaH, 'LineWidth', 2);
plot(time, -(JFa - JFaH), 'LineWidth', 2)
hold off
grid on
xlabel('Days')
ylabel('\Delta Q');
legend('-J_F^z', '-J_{F-Horiz}^z', '-J_{F-Vert}^z');
set(gca, 'FontSize', fs);
set(gcf, 'Color', 'w');

subplot(2,1,2)
plot(time, -JFa, 'LineWidth', 2);
hold on
plot(time, -JBa, 'LineWidth', 2);
plot(time, -JFaN, 'LineWidth', 2);
plot(time, -JBaN, 'LineWidth', 2)
% plot(time, -JFaH, 'LineWidth', 2);

hold off
grid on
xlabel('Days')
ylabel('\Delta Q');
legend('-J_F^z','-J_B^z', '-J_{F-Num}^z', '-J_{B-Num}^z');
set(gca, 'FontSize', fs);
set(gcf, 'Color', 'w');