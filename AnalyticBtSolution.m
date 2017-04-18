%% ANALYTIC D_T solution

ebf = 1;
sigma = .8;

t=linspace(.1, 10, 100).';

ef = double(erfc(1,sym(0)));
% irfc = 
B_t = sqrt(2*t*sigma)*2.*ebf.*ef;
dBdt = sqrt(2*sigma)./sqrt(t).*ebf.*ef;
dBdt = sqrt(2*sigma)./sqrt(pi*t).*ebf;
dBdt = 1./(sqrt(pi*t./2)).*ebf; % Time normalized by t = t'/f;
dBdt = 1./(pi*sqrt(t)).*ebf; % Normalized by t = 2*pi*t'/f ie. inertial period.
dBdtd = gradient(B_t, t(2)-t(1));
D = -dBdt;
% subplot(1,2,1)
% plot(t, B_t);
% subplot(1,2,2)

figure
subtightplot(1,1, 1, [1 1], .2, .25)
plot(t, D, 'LineWidth', 2, 'Color', 'k')
hold on;
% plot(t, -dBdtd, '--');
hold off
% hold on;
% plot(t, dBdt, '--');
% hold off
grid on
set(gca, 'FontSize', 16)
xlabel('$\frac{t}{(2 \pi/f)}$', 'FontSize', 20)
ylabel('$\frac{J_{D_t}}{(f\;EBF_g/\hat{\delta_t})}$', 'Rotation', 0, 'FontSize', 20)

set(gca, 'ylim', [-1 0]);
set(gcf, 'Color', 'w', 'Position', [ 670         499        1001         472]);
%% TEST the solution to the equation
sigma = 0.5;

z = linspace(-10, 0, 102);
bf = NaN(length(t), length(z));
bzz = bf;
for i=1:length(t)
    disp(i);
    zn = z./sqrt(2*t(i)./sigma);
    
    efs = erfc(1,sym(-zn));
    ef = double(efs);
%     efz = diff(efs, sym(zn));
%     efzz = diff(efs, sym(zn)), 
    bf(i,:) = sqrt(2*t(i)./sigma).*ebf.*ef;
%     bz(i,:) = diff(
    bzz(i,:) = 4*del2(bf(i,:), z(2)-z(1))./(2*sigma);
end




[bz, bt] = gradient(bf,  z(2)-z(1), t(2)-t(1));
% [bzz, ~] = gradient(bz,  z(2)-z(1), 1e10);
% for i=1:length(
% bzz = 4*del2(bf, z(2)-z(1));
%%

d  = length(z)-0;
subplot(3,1,1)
plot(bt(:,d));
hold on
plot(bzz(:,d), '--')
hold off;

subplot(3,1,2)
plot(bz(:,end)./ebf);
set(gca, 'ylim', [0.5 1.5])

subplot(3,1,3)
plot(bt(10,:), z)
hold on
plot(bzz(10,:), z, '--');
hold off
%%
pcolor(t, z, bf); shading interp; colorbar;

