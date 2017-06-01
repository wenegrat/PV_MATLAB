%%
nl = length(conts);
haxnew = NaN(nl-1, 3);
haxold = cptcmap('GMT_haxby');
for i=1:3
    haxnew(:,i) = interp1(linspace(conts(1), conts(end),32), haxold(:,i), linspace(conts(1), conts(end),nl-1));
end
%%
gradbs = logspace(-8, -6, 100);
gradbs = linspace(1e-8, 1e-6, 1000);
hs = linspace(20, 500, 500);
clear vals
f = 1e-4;
for i=1:length(hs)
vals(i,:) = hs(i)^2.*0.05.*(gradbs).^2./(f).*1035*3994./(9.81*2e-4);
end
conts = linspace(0, 500, 20);
conts = 0:50:1000;

figure
subtightplot(1,1,1, [0 0], .2, .2);
[c, h] =contourf(hs, gradbs./f^2, vals.', conts);
% pcolor(hs, gradbs./f^2, vals.'); shading interp


set(gca, 'YScale', 'log')
xt = get(gca, 'XTick');
% hold on
% plot(xt, ones(size(xt)).*(2*(2*2*pi.*sind(38)./86400))^2)
% plot(xt, ones(size(xt)).*(4*(2*2*pi.*sind(38)./86400))^2)
% plot(xt, ones(size(xt)).*(6*(2*2*pi.*sind(38)./86400))^2)
% hold off
% clabel(c)
grid on
xlabel('$H$');
set(gca, 'FontSize', 20)

ylabel('$\frac{|\nabla_h b|}{f^2}$', 'Rotation', 0, 'FontSize', 26);
colormap(haxnew)
title('Equivalent surface heat flux');
cb = colorbar;
set(get(cb, 'yLabel'), 'String',  '$\mathrm{W m^{-2}}$', 'Interpreter', 'Latex', 'Rotation', 0);
set(cb, 'Ticks', 0:200:1000, 'TickLabels', {'0', '200', '400', '600', '800', '$>$1000'});
set(gcf, 'Color', 'w', 'Position', [          296         435        1092         477]);

latmixf = 2*2*pi./(86400).*sind(38)
obsbuoy = 5e-7
obsbuoy = 1e-8;
pos = obsbuoy./latmixf^2
h = 250;
hold on
% plot(h, pos, 'x');
hold off
%%
clear g
f = 1e-4;
M2 = [(1*f)^2 (2*f)^2 (4*f)^2 (6*f)^2];
Q = [25 100 200];
% M2 = linspace(f^2, (6*f)^2, 90);
% Q = linspace(25, 200, 100);
h = 150;
disp(' ')
disp(' ')
disp('=======================')
for i=1:length(M2)
    for j=1:length(Q)
        disp('--------------------');
        disp(['M^2_o = ', num2str(M2(i))]);
        disp(['Q_o  = ', num2str(Q(j))]);
        QtoB = (9.81*2e-4./(1035*3994));
        Qe = h^2*0.05*M2(i)^2./(f) * 1./QtoB;
        disp(['Q_{Eff}= ', num2str(Qe)]);
        Bo = QtoB.*Q(j);
        gamma = 0.05*h^2*M2(i).^2./(f*Bo);
        disp(['Gamma  = ', num2str(gamma)]);
        Ri = (64*f)^2*f^2./(M2(i).^2);
       
        taustone = sqrt(54/5)*sqrt(1+Ri)./f; % See FK08 eq 3
        disp(['Ri = ', num2str(Ri)]);
        disp(['ts = ', num2str(taustone./(86400))]);
        
       
        g(i,j) = gamma;
    end
end
%%
contourf(Q, M2, log10(g),30); xlabel('Q'); ylabel('M2');
colorbar