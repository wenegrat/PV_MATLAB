%% Thomas and Ferrari Equation 6
statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc'; extrafile = 'extra.nc';
kppfile = 'kppdiags.nc';

f0 = 1e-4;

% hl = 36;
Z = ncread(statefile, 'Z');
Zl = ncread(statefile, 'Zl');
nx = 150; ny =200;
nx = 3;
% metric = permute(repmat(Z,  [1, nx, ny, nt]), [2 3 1 4]);

dx = 500;
dy = 500;
ts = 7200;
st = 1;
nt = 539;
nt = 719;
metric = permute(repmat(Z,  [1, nx, ny, nt]), [2 3 1 4]);

slice = {0, 0, 0, [st st+nt-1]};
b = GetVar(statefile, diagfile, {'b', '(1)'},slice);
%%
bx = DPeriodic(b, dx, 'x');
by = DPeriodic(b, dy, 'y');
bz = Drv(metric, b, 'z');
ZETA = GetVar(statefile, extrafile, {'momVort3', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);
%%
% bzbar = squeeze(nanmean(nanmean((b(:,:,1,:) - b(:,:,hl,:))./(Ho))));
bbar = squeeze(nanmean(nanmean(b)));
% m2bar = squeeze(nanmean(nanmean( abs(bx(:,:,2,:)+1i.*by(:,:,2,:)))));
% iribar = m2bar.^2./(f0.^2.*bzbar);
%%
H = GetVar(statefile, etanfile, {'KPPhbl', '(1)'}, {slice{1},slice{2}, [1 1], slice{4}});

[bt, ~] = gradient(bbar, ts);

zetab = ZETA.*b;
zba = squeeze(nanmean(nanmean(zetab)));
% zba = squeeze(zba(1,:) - zba(hl,:));
[zbat, ~] = gradient(zba, ts);

% zbat = -1./(f0*Ho).*zbat;


wq = W.*outputFull.Q(:,:,:,st:(st+nt-1));
wqa = squeeze(nanmean(nanmean(wq)));
% wqa = squeeze(wqa(2,:) - wqa(hl,:));
% wqa = -1./(f0*Ho).*wqa;


Jf = outputFull.JFz(:,:,:,st:(st+nt-1));
Jfa = squeeze(nanmean(nanmean(Jf)));
% Jfa = squeeze(Jfa(2,:) - Jfa(hl,:));
% Jfa = -1./(f0.*Ho).*Jfa;

Jb = outputFull.JBz(:,:,:,st:(st+nt-1));
Jba = squeeze(nanmean(nanmean(Jb)));
% Jba = squeeze(Jba(2,:) - Jba(hl,:));
% Jba = -1./(f0.*Ho).*Jba;

% test= outputFull.JBz(:,:,:,st:(st+nt-1));
% test = test.*f0./(f0+ZETA);
% test = squeeze(nanmean(nanmean(test)));
% test = squeeze(test(2,:) - test(hl,:));
% test = -1./(f0*Ho).*test;
% %%
% % SCALE = squeeze(nanmean(nanmean(repmat(H, [1 1 nz 1]).*abs(bx+1i.*by).^2)));
% % SCALE = squeeze(nanmean(nanmean(abs(bx+1i.*by).^2)))./f0;
% % 
% % SCALE = squeeze(nanmean(SCALE(1:hl,:)));
% dmudz0 = -104./(21.*output.H);
% dmudzHo = -8*(output.H - 2*Ho).*(20*Ho.^2 - 20.*output.H.*Ho + 13.*output.H.^2)./(21*output.H.^4);
% dmudzHo(output.H<Ho) = 0;
% dmudzHo = 0; %xx- not formally correct;
% % FKSCALE = 0.6.*f0.*output.MagB;
% FKSCALE = -0.06.*output.H.^2.*output.MagB.^2.*(dmudz0-dmudzHo)./(f0.*Ho);
% % FKSCALE = -0.06.*output.H.^2.*output.MagB.*(dmudz0-dmudzHo)./(f0.*Ho).*squeeze(nanmean(nanmean(nanmean(abs(bx(:,:,1:hl,:) + 1i.*by(:,:,1:hl,:))))));
% % FKSCALE = Nan(length(FKSCALE),1);
% % FKSCALE(Ho>output.H) = 0.3.*f0.*output.MagB(Ho>output.H).*output.H(Ho>output.H)./Ho;
% SURFSCALE = -1./(f0.*Ho).*output.dJbsa.*1.2./(200*180*500*500);
%%
zlu = 2;
zld = 22;16;
Ho = abs(Z(zld));

figure
subplot(2,1,1)

plot(outputFull.time, (bt(zlu,:) - bt(zld,:))./Ho, 'LineWidth', 3);
hold on
plot(outputFull.time,-(zbat(zlu,:)-zbat(zld,:))./(f0.*Ho), 'LineWidth', 2);
plot(outputFull.time,-(wqa(zlu,:)-wqa(zld,:))./(f0.*Ho), 'LineWidth', 2);
plot(outputFull.time,-(Jba(zlu,:)-Jba(zld,:)+Jfa(zlu,:) - Jfa(zld,:))./(f0.*Ho), '-', 'LineWidth', 2);
hold off
set(gca, 'FontSize', 16, 'xlim', [0 45], 'ylim', [-3 4].*1e-12)

% Label dN^2/dt
[x1 y1] = ds2nfu(8, 2.1e-12);
[x2 y2] = ds2nfu(10, 1.6e-12);
annotation('textarrow', [x1 x2], [y1 y2], 'String', '$\frac{\partial \overline{N^2}}{\partial t}$',...
    'Interpreter', 'Latex', 'FontSize', 24, 'Color', [  0    0.4470    0.7410])
% Label FRONT
[x1 y1] = ds2nfu(15, -1.4e-12);
[x2 y2] = ds2nfu(16, -0.65e-12);
annotation('textarrow', [x1 x2], [y1 y2], 'String', '$FRONT$',...
    'Interpreter', 'Latex', 'FontSize', 20, 'Color', [0.8500    0.3250    0.0980])
% Label WQ
[x1 y1] = ds2nfu(8, -0.8e-12);
[x2 y2] = ds2nfu(9, 0e-12);
annotation('textarrow', [x1 x2], [y1 y2], 'String', '$ADV$',...
    'Interpreter', 'Latex', 'FontSize', 20, 'Color', [0.9290    0.6940    0.1250])
% Label NC
[x1 y1] = ds2nfu(23, 2.1e-12);
[x2 y2] = ds2nfu(21.5, 1.7e-12);
annotation('textarrow', [x1 x2], [y1 y2], 'String', '$JNC$',...
    'Interpreter', 'Latex', 'FontSize', 20, 'Color', [0.4940    0.1840    0.5560])
% legend('$\frac{\partial \overline{N^2}}{\partial t}$', '$-\frac{1}{fH}\frac{\partial <\zeta b>}{\partial t}$', '$-\frac{1}{fH}<wq>$',...
%     '$-\frac{1}{fH}<J_F^z>$', '$-\frac{1}{fH}<J_D^z>$', '$-\frac{1}{fH}(<J_F^z>+<J_D^z>)$')
% title(['z = -(' num2str(abs(Z(zlu)),2),' - ',num2str(abs(Z(zld)),2), ' m)']);
title(['$z_t = -$' num2str(abs(Z(zlu)),2),' m,       $z_b= -$',num2str(abs(Z(zld)),2), ' m']);

grid on
xlabel('Days');
ylabel('$s^{-3}$');

zlu = 10;2;
zld = 17;36;
Ho = abs(Z(zld));

subplot(2,1,2)

plot(outputFull.time, (bt(zlu,:) - bt(zld,:))./Ho, 'LineWidth', 3);
hold on
plot(outputFull.time,-(zbat(zlu,:)-zbat(zld,:))./(f0.*Ho), 'LineWidth', 2);
plot(outputFull.time,-(wqa(zlu,:)-wqa(zld,:))./(f0.*Ho), 'LineWidth', 2);
plot(outputFull.time,-(Jba(zlu,:)-Jba(zld,:)+Jfa(zlu,:) - Jfa(zld,:))./(f0.*Ho), '-', 'LineWidth', 2);
hold off
set(gca, 'FontSize', 16, 'xlim', [0 45], 'ylim', [-3 4].*1e-12)

% Label dN^2/dt
[x1 y1] = ds2nfu(8, 2.1e-12);
[x2 y2] = ds2nfu(9, 1.8e-12);
annotation('textarrow', [x1 x2], [y1 y2], 'String', '$\frac{\partial \overline{N^2}}{\partial t}$',...
    'Interpreter', 'Latex', 'FontSize', 24, 'Color', [  0    0.4470    0.7410])
% Label FRONT
[x1 y1] = ds2nfu(17, -1.3e-12);
[x2 y2] = ds2nfu(16.25, -0.5e-12);
annotation('textarrow', [x1 x2], [y1 y2], 'String', '$FRONT$',...
    'Interpreter', 'Latex', 'FontSize', 20, 'Color', [0.8500    0.3250    0.0980])
% Label WQ
[x1 y1] = ds2nfu(24, 2.1e-12);
[x2 y2] = ds2nfu(23, 1.5e-12);
annotation('textarrow', [x1 x2], [y1 y2], 'String', '$ADV$',...
    'Interpreter', 'Latex', 'FontSize', 20, 'Color', [0.9290    0.6940    0.1250])
% Label NC
[x1 y1] = ds2nfu(7, -1.1e-12);
[x2 y2] = ds2nfu(7.5, -0.35e-12);
annotation('textarrow', [x1 x2], [y1 y2], 'String', '$JNC$',...
    'Interpreter', 'Latex', 'FontSize', 20, 'Color', [0.4940    0.1840    0.5560])
% legend('$\frac{\partial \overline{N^2}}{\partial t}$', '$-\frac{1}{fH}\frac{\partial <\zeta b>}{\partial t}$', '$-\frac{1}{fH}<wq>$',...
%     '$-\frac{1}{fH}<J_F^z>$', '$-\frac{1}{fH}<J_D^z>$', '$-\frac{1}{fH}(<J_F^z>+<J_D^z>)$')
title(['$z_t = -$' num2str(abs(Z(zlu)),2),' m,       $z_b= -$',num2str(abs(Z(zld)),3), ' m']);

grid on
xlabel('Days');
ylabel('$s^{-3}$');

set(gcf, 'Color', 'w', 'Position', [680          92        1033         830])