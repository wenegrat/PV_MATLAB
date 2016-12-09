gradbs = logspace(-8, -6, 100);
hs = linspace(20, 300, 500);
clear vals
for i=1:length(hs)
vals(i,:) = hs(i)^2.*0.04.*(gradbs).^2./(2*2*pi.*sind(38)./86400).*1035*3994./(9.81*2e-4);
end
conts = linspace(0, 500, 20);
[c, h] =contour(hs, gradbs, vals.', conts);
set(gca, 'YScale', 'log')
clabel(c)
