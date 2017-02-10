

%% First just look at how closed time integrated PV budget is.
figure
subplot(2,1,1)
plot(output.Qa); 
hold on
plot(-output.Jfa);
plot(-output.Jba);
plot(-(output.Jfa + output.Jba));
hold off

cr = corr(output.Qa, -(output.Jfa+output.Jba));

r = regress(output.Qa, -(output.Jfa+output.Jba));

title(['Corr: ', num2str(cr), '   Regress Coeff: ', num2str(r)]);

subplot(2,1,2)
plot(output.Qt); 
hold on
plot(-output.dJf);
plot(-output.dJb);
plot(-(output.dJf + output.dJb));
hold off

cr = corr(output.Qa, -(output.Jfa+output.Jba));

r = regress(output.Qa, -(output.Jfa+output.Jba));

title(['Corr: ', num2str(cr), '   Regress Coeff: ', num2str(r)]);

%% Scatter Plots
figure
% VERIFY BASIC PV BUDGET
subplot(1,2,1)
scatter(output.Qt, -(output.dJf+output.dJb));
grid on
cr = corr(output.Qt, -(output.dJf + output.dJb));
r = regress(output.Qt, -(output.dJf+output.dJb));

title(['Corr: ', num2str(cr), '   Regress Coeff: ', num2str(r)]);
xlabel('dQt/dt');
onetoone

subplot(1,2,2)
scatter(output.Qa, -(output.Jfa+output.Jba));
grid on
cr = corr(output.Qa, -(output.Jfa + output.Jba));
r = regress(output.Qa, -(output.Jfa+output.Jba));

title(['Corr: ', num2str(cr), '   Regress Coeff: ', num2str(r)]);
xlabel('Delta Q');
onetoone

