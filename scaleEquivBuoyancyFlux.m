gradbs = logspace(-8, -6, 100);
hs = linspace(20, 300, 500);
clear vals
for i=1:length(hs)
vals(i,:) = hs(i)^2.*0.04.*(gradbs).^2./(2*2*pi.*sind(38)./86400).*1035*3994./(9.81*2e-4);
end
conts = linspace(0, 500, 50);
[c, h] =contour(hs, gradbs, vals.', conts);
set(gca, 'YScale', 'log')
xt = get(gca, 'XTick');
hold on
plot(xt, ones(size(xt)).*(2*(2*2*pi.*sind(38)./86400))^2)
plot(xt, ones(size(xt)).*(4*(2*2*pi.*sind(38)./86400))^2)
plot(xt, ones(size(xt)).*(6*(2*2*pi.*sind(38)./86400))^2)
hold off
clabel(c)


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
        disp(['M2 = ', num2str(M2(i))]);
        disp(['Q  = ', num2str(Q(j))]);
        QtoB = (9.81*2e-4./(1035*3994));
        Qe = h^2*0.04*M2(i)^2./(f) * 1./QtoB;
        disp(['Q_E= ', num2str(Qe)]);
        Bo = QtoB.*Q(j);
        gamma = 0.04*h^2*M2(i).^2./(f*Bo);
        disp(['G  = ', num2str(gamma)]);
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