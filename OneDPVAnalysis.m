statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc'; extrafile = 'extra.nc';
kppfile = 'kppdiags.nc';

TtoB = 9.81.*2e-4;
f0 = 1e-4;

Q0 = squeeze(ncread(etanfile, 'TFLUX'));
QNET = squeeze(ncread(etanfile, 'oceQnet'));
QSW = squeeze(ncread(etanfile, 'oceQsw'));
X = ncread(statefile, 'X');
Y = ncread(statefile, 'Y');
Z = squeeze(ncread(statefile, 'Z'));
Zl = ncread(statefile, 'Zl');
T = ncread(diagfile, 'T');
dh = diff([Zl; -300]);
dh = abs(dh);

nx = length(X); ny = length(Y); nz = length(Z);
% dx = X(2)-X(1)
% dy = Y(2)-Y(1)
dz = Z(1)-Z(2) %surface only, XX-should track this through the code to ensure correct.
ts = T(2)-T(1)
slice = {0, 0, 0, 0};

b = TtoB.*squeeze(ncread(diagfile, 'THETA'));
[nz nt] = size(b);
bz = (b(1:end-1,:) - b(2:end,:))./repmat((Z(1:end-1) - Z(2:end)), [1 nt]);
% bz =  Drv(Zl, b, 'z');
% bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);

q = f0.*bz;

Qa = nansum(q.*repmat(dh(1:end-1), [1 nt]));
Qa = Qa - Qa(1);
JD = gradient(Qa, ts);
cd
h = squeeze(ncread(etanfile, 'KPPhbl'));
Bo = 9.81*2e-4*squeeze(Q0(:))./(1035*3994);
Bsw = 9.81*2e-4.*squeeze(QSW(:))./(1035*3994);
% h(h<6) = 6;
h(h<3) = 3;
BLW = Bo-Bsw;
r = .58;
d1 = 0.35;
d2 = 23;
hs = h;
Bsw1 = r*Bsw.*(1-exp(-3./d1));
Bsw2 = (1-r).*Bsw.*(1-exp(-3./d2));
% Bsw2 = ;
% Bsw1 = r.*Bsw;
JS = f0.*Bo./h;
JSW = JS;
JSW = f0./hs.*(BLW+Bsw1+Bsw2);
% JSW = f0.*(Bsw1./3+ Bsw2./3 + BLW./h) ;
% JSW = 0.057.*f0./h.^2.*(BLW.*h./2 + r*Bsw.*(3+h./2) + (1-r).*Bsw.*(23+h/2));
% JSW = f0.*((BLW)./3 +  r.*Bsw./3 + (1-r).*Bsw*(1-exp(-3./d2))./3);
try
tx = ncread(etanfile, 'oceTAUX');
tx = squeeze(tx(1,:,:,:));
JFWe = 0.7*sqrt(tx./1035).^3.*1e-4./(h.^2);
catch
    disp('here');
    JFWe = 0;
    tx = 0;
end
ht = gradient(h);
ustar = sqrt(tx./1035);
hobu = ustar.^3./(0.4.*Bo);
he = 0.7.*ustar./f0;
wstar = (-Bo.*h).^(1/3);
% JFWe = 0.7.*( ustar.^3 +(-Bo.*h)).*f0./(h.^2);
% JFWe( (ustar.^3 + (-Bo.*h))<0) = 0;
% JSW = f0.*(.1.*h.*ustar.*( (r*Bsw)./3.^2 + (1-r).*Bsw./d2.^2)  + r*Bsw./3 + (1-r).*Bsw./d2);
% JSW = h.*ustar.*(Bo);
% JFWe(ht<0) = 0;
% JSW = JSW+ JFWe;
% JSW(BLW+Bsw1+Bsw2 > 0) = JS(BLW+Bsw1+Bsw2>0);
% JSWT = f0./h.*(BLW+Bsw1+Bsw2);
% JSW((BLW+Bsw1+Bsw2)>0) = JSWT((BLW+Bsw1+Bsw2)>0);
% JSWT = f0./h.*(BLW+Bsw1+Bsw2);
% JSW((BLW+Bsw1+Bsw2)<0) = JSWT((BLW+Bsw1+Bsw2)<0);

% JSW = f0./3.*(BLW+Bsw1+Bsw2);% + (1-0.58).*Bsw./h.*(1-exp(-3/23)));

% JSW = JSW - f0./h.*((
%%
xl = 1:1200;
figure
subplot(3,2,1:2)

plot(JD(xl), 'LineWidth', 2);
hold on
plot(JSW(xl));
plot(JS(xl), '--');
hold off
legend('J_{D_{NUM}}', 'J_{S-GAIN}', 'J_{S-LOSS}');
grid on

subplot(3,2,3:4)
plot(-h(xl))
hold on
plot(-hobu(xl));
% plot(-he(xl));
hold off
grid on
ylabel('h_{KPP}');
set(gca, 'ylim', [-300 0]);

subplot(3,2,5)
posl =50:600;
scatter(JD(posl), JS(posl));
onetoone
rc = regress(JD(posl).', JS(posl));
cr = corr(JD(posl).', JS(posl));

hold on
hold off
title(['Time Period: ', num2str(posl(1)), '-', num2str(posl(end)),'  Regress Coeff: ', num2str(rc,2), '   Corr: ', num2str(cr,2)]);
grid on

subplot(3,2,6)
posl =850:1000;
scatter(JD(posl), JSW(posl));
onetoone
rc = regress(JD(posl).', JSW(posl));
cr = corr(JD(posl).', JSW(posl));

hold on
hold off
title(['Time Period: ', num2str(posl(1)), '-', num2str(posl(end)),'  Regress Coeff: ', num2str(rc,2), '   Corr: ', num2str(cr,2)]);
grid on

set(gcf, 'Color', 'w', 'Position', [675         196        1049         778]);

%%
% subplot(3,1,3)
figure
plot(BLW(xl));
 hold on
 plot(Bsw(xl));
 plot(BLW(xl) + Bsw(xl), 'LineWidth', 2);
 plot(Bsw1(xl), '--');
 plot(Bsw2(xl),'--');
 plot(BLW(xl) + Bsw1(xl) );
 hold off
 grid on
 legend('B_{LW}', 'B_{SW}', 'Bo');
% plot(Bsw1(xl));
% plot(Bsw2(xl));
% plot(BLW(xl)+Bsw1(xl)+Bsw2(xl));
% plot(Bo(xl), '--');
% hold off
% plot(BL
% subplot(3,1,3)
% plot(Q0)
% hold on
% plot(Q0-QSW);
% plot(QSW);
% plot(Q0 - QSW + 0.58*QSW); % Fluxes at surface bin (LW - r1*SW)
% hold off
% grid on
%%
r = -b./(g./1035);
mlmask = r > repmat(r(1,:), [50 1]) + 0.005;
mlmask = double(mlmask);
mlmask(mlmask==0) = NaN;
zf = repmat(Z , [1 nt]);
mldepth =squeeze( max(zf.*mlmask, [], 1));
%%
xl = 100:1000;
zl = 1:40;
figure
pcolor(xl, squeeze(Z(zl)), squeeze(bz(zl,xl)));shading interp
hold on
plot(xl, -squeeze(h(xl))); 
plot(xl, mldepth(xl));
hold off
%%
figure
subplot(2,1,1)
posl = 10:525;
scatter(JD(posl), JS(posl));
onetoone
rc = regress(JD(posl).', JS(posl));
hold on
% scatter(JD(posl), -rc.*JS(posl), 'x');
hold off
title(num2str(rc))

subplot(2,1,2)
posl = 850:950;
scatter(JD(posl), JSW(posl));
onetoone
title(num2str(regress(JD(posl).', JSW(posl))))

%%
posl =600:800;
scatter(-1./h(posl),(JD(posl) -JSW(posl).')./(f0.*Bo(posl).'))
onetoone
regress(-1./h(posl), ((JD(posl) - JSW(posl).')./(-f0.*Bo(posl).')).')
grid on
% axis equal