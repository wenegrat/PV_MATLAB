%% ITERATIVELY LOAD OUTPUT FILES
% output.Qa = Qa;
% output.Qt = Qt;
% output.dJf = dJFdt;
% output.dJb = dJBdt;
% output.Jfa = JFa;
% output.Jba = JBa;
% output.dJfa_t = Jftota;
% output.dJba_t = Jbtotpa;
% output.dJbsa = Jbsurfa;
% output.dJbea = Jbeddya;
% output.Tsurf = squeeze(THETA(:,:,1,:));
% output.Q = Q0(1,1,1,1);
% output.time = 
inputs = output;

%% MAKE FIGURES
figure
% VERIFY BASIC PV BUDGET
subplot(1,2,1)
scatter(inputs.Qt, -(inputs.dJf+inputs.dJb));
grid on
cr = corr(inputs.Qt, -(inputs.dJf + inputs.dJb));
title(num2str(cr));
xlabel('dQt/dt');

subplot(1,2,2)
scatter(inputs.Qa, -(inputs.Jfa+inputs.Jba));
grid on
cr = corr(inputs.Qa, -(inputs.Jfa + inputs.Jba));
title(num2str(cr));
xlabel('Delta Q');

figure
% Verify Theory Scalings
subplot(1,2,1)
scatter(-inputs.dJf, -inputs.dJfa_t);
grid on
cr = corr(-inputs.dJf, -(inputs.dJfa_t));
title(num2str(cr));
subplot(1,2,2)
scatter(-inputs.dJb, -inputs.dJba_t);
grid on
cr = corr(-inputs.dJb(3:end), -(inputs.dJba_t(3:end)));
title(num2str(cr));
%%
figure
subplot(2,1,1)
% indj = 30; indk = 76;
% plotyy(time, abs(squeeze(JFz(indj, indk, 2,:))), time, squeeze(JfTtot(indj, indk,:)));
% [ax, h1, h2] = plotyy(time, -inputs.dJf , time, -inputs.dJfa_t);
% set(h1, 'LineWidth', 2);
% set(h2, 'LineWidth', 2);
plot(time, -inputs.dJf, 'LineWidth', 2);
hold on
plot(time, -inputs.dJfa_t, 'LineWidth', 2);
hold off
grid on
% title(num2str(corr(dJFdt, Jftota)))
cr = num2str(corr(inputs.dJf, inputs.dJfa_t));
Bf = (regress(inputs.dJf, inputs.dJfa_t));
title(['Corr: ', cr, '  Regress Coeff: ', num2str(Bf)])

subplot(2,1,2)
jvar = inputs.dJb;
% [ax, h1, h2]=plotyy(time, -jvar, time, -inputs.dJba_t);
% set(h1, 'LineWidth', 2);
% set(h2, 'LineWidth', 2);
plot(time, -jvar, 'LineWidth', 2);
grid on
hold on
plot(time, -inputs.dJba_t, 'LineWidth', 2);
% plot(time, -2*Jbtotpa);
% plot(time, Jbtotpa, '--');
hold off
cr = num2str(corr(jvar(1:end), inputs.dJba_t(1:end)));
Bb = (regress(jvar, inputs.dJba_t));
title(['Corr: ', cr, '  Regress Coeff: ', num2str(Bb)])

% subplot(3,1,3)
% % plotyy(time, abs(dJFdt./dJBdt), time, abs(Jftota./Jbtota));
% plot(time, (-dJFdt+dJBdt), 'LineWidth',2);
% hold on
% plot(time,(-Jftota*Bf+Jbtotpa*Bb), 'LineWidth', 2);
% hold off
% grid on

%%
tind = 100;
isoT = [0 100];
Qa = inputs.Qa;
JFa = inputs.Jfa;
JBa = inputs.Jba;
Qt = inputs.Qt;
dJFdt = inputs.dJf;
dJBdt = inputs.dJb;
titleString = ['Full Volume           Surface B_0: ', num2str(squeeze(Q0(1,1,1)))];
QBudgetPlot;
dQdtPlot;