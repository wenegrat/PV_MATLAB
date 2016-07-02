%%
dx = 1000; dy = 1000; dz = 2.5;
TtoB = 9.81.*2e-4;
% TtoB = 1;
% Following JAS article
bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);
THETA = GetVar(statefile, diagfile, {'THETA', '(1)'},slice);
b = TtoB.*THETA;
%%
% Equation 1b
% Is my vorticity incompressible?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OMEGAX = GetVar(statefile, diagfile, {'WVEL', 'VVEL', 'Dy(1) - Dz(2)'}, slice);
% OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'WVEL', 'Dz(1) - Dx(2)'}, slice);
OMEGAX = GetVar(statefile, diagfile, {'VVEL', ' - Dz(1)'}, slice);
OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);
OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL', '(1) + Dx(2) - Dy(3)'}, slice);
ztmp = ncread(statefile, 'Z');
[nx, ny, nz, nt] = size(OMEGAX);
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);
div = Drv(dx, b.*OMEGAX, 'x') + Drv(dy, b.*OMEGAY, 'y') + Drv(metric, b.*OMEGAZ, 'z');
Q = bx.*OMEGAX + by.*OMEGAY + bz.*OMEGAZ;

res1b = (Q - div)./Q;
% Find residual is not always small locally, but seems small in a volume
% integrated budget. Normalized mean is O(-3e-3). Although some instances
% of larger values...
figure('Name', 'Eq. 1b residual')
plot(squeeze(nanmean(nanmean(nanmean(res1b)))))

%%
% Check 2A, note that there is an additional term (Q/rho drho/dt)
bt = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND', '(1)/86400'}, slice);
Term2 = 1./9.81.*Q.*bt; %If this value is O(1/ts) then it could matter
Term2a = averageQForcings(Term2, mask, gridvol, ts);
figure('Name' , 'Eq. 2a extra term')
plot(Term2a./Qa)
% Note I'm not sure if this should be in terms of buoyancy or what, units
% are wrong here...
% scatter(Qa+Fricst+Diast, Term2a)
grid on
%%
% Ensure that baroclinic term disappears (Eq. 4) when dotted with buoyancy
% gradient.

% baroi = by.*b*1035 - bz.*dpdy;
% baroj = bz.*dpdx - bx.*b*1035;
% barok = bx.*dpdy - by.*dpdx;
% 
% dotted = 1/(1035^2).*(bx.*by.*b*1035 -bx.*bz.*dpdy + by.*bz.*dpdx -by.*bx.*b*1035 + bz.*bx.*dpdy -bz.*by.*dpdx);
%Ok, cancels regardless of operators.
%%
% Equation 6a
% Can I balance my budget like this?
% Which Q is balanced here (direct or normal)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frictional Terms
dpdx = GetVar(statefile,diagfile, {'Um_dPHdx','(1)'},slice);
dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.8*Dx(1)'},sliceEta));
dpdx = permute(repmat((dEtadx), [1, 1, 1, nz]), [1 2 4 3]) + dpdx;
dpdy = GetVar(statefile,diagfile, {'Vm_dPHdy','(1)'},slice);
dEtady = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.8*Dy(1)'},sliceEta));  
dpdy = permute(repmat((dEtady), [1, 1, 1, nz]), [1 2 4 3]) + dpdy;

Fx = GetVar(statefile, diagfile, {'TOTUTEND','Um_Advec', ' (1)/86400 - (2)'},slice);
Fx = Fx - dpdx;
Fy = GetVar(statefile, diagfile, {'TOTVTEND', 'Vm_Advec', '(1)/86400 - (2)'}, slice);
Fy = Fy - dpdy;

Fyu = Fy;
Fyu(:,end-1:end,:,:) = 0; %Ad hoc removal of gradient at wall.
dFxdz = Drv(metric, Fx, 'z');
dFxdy = Drv(dy, Fx, 'y');
dFydz = Drv(metric, Fyu, 'z');
dFydx = Drv(dx, Fyu, 'y');
byu = by;
byu(:,end-1:end,:,:) = 0;
bzu = bz;
bzu(:,:,end,:) =0;

FrictTerm = -bx.*dFydz + byu.*dFxdz + bzu.*(dFydx - dFxdy);
%Diabatic Terms
D =  TtoB*GetVar(statefile, diagfile, {'TOTTTEND', 'UDIAG1', '(1)/86400 - (2)'}, slice);
Dx = Drv(dx, D, 'x');
Dy = Drv(dy, D, 'y');
Dz = Drv(metric, D, 'z');

DiaTerm = OMEGAX.*Dx + OMEGAY.*Dy + OMEGAZ.*Dz;

FT = averageQForcings(FrictTerm, mask, gridvol, ts);
DT = averageQForcings(DiaTerm, mask, gridvol, ts);

%Advective Terms
Q = bx.*OMEGAX + byu.*OMEGAY + bzu.*OMEGAZ;
Qx = Drv(dx, Q, 'x');
Qy = Drv(dy, Q, 'y');
Qz = Drv(metric, Q, 'z');
U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

AdvTerm =  -U.*Qx - V.*Qy- W.*Qz;
AT = averageQForcings(AdvTerm, mask, gridvol, ts);
%%
figure('Name', 'Eq 6a balance');
plot(Qa);
hold on
plot(Qda, '-.');
plot(FT);
plot(DT);
plot(AT);
plot(FT+DT,'--');
hold off
legend('Q', 'Qdir', 'Fric', 'Dia', 'Adv', 'Sum');
%Note that Frictional term seems very strongly affected by the wall
%condition. Maybe need to go back and look at how this is being calculated?
%Or consider doing doubly periodic domain? This would be consistent with
%Thomas and Ferrari 2008. Can I explain why the 'direct' Q calc would
%include this effect?

%%
% Look at Eq. 12, cancellation between WAE and MAE?
%%
% Eq 13, not that w~=0 at surface.
ADs = squeeze(nansum(nansum(JAz(:,yl:end-yl,2,:).*mask(:,yl:end-yl,2,:))) - nansum(nansum(JAz(:,yl:end-yl,end-zl+1,:).*mask(:,yl:end-yl,end-zl+1,:)))).*dx.*dy;
ADst = cumtrapz(ADs).*ts./vol;
%O(1e-13) too small.

%%
% They recommend using Eq. 15 as a  first diagnostic. try that here. Wait
% though, this is supposed to be outside the BL?
%FRICTION TERMS
dFudz = Drv(metric, Fx, 'z');
dFvdz = Drv(metric, Fy, 'z');
dFudy = Drv(dy, Fx, 'y');
dFvdx = Drv(dx, Fy, 'x');
WF = bx.*(-dFudz) + by.*(dFvdz) + bz.*(dFvdx - dFudy);
WFa = averageQForcings(WF, mask, gridvol, ts);

iterm = bx.*(Drv(dy, OMEGAX.*V - OMEGAY.*U, 'y') - Drv(metric, OMEGAZ.*U - OMEGAX.*W,'z'));
jterm = by.*(Drv(metric,OMEGAY.*W - OMEGAZ.*V,'z') - Drv(dx, OMEGAX.*V - OMEGAY.*U, 'x'));
kterm = bz.*(Drv(dx, OMEGAZ.*U - OMEGAX.*W, 'x') - Drv(dy, OMEGAY.*W - OMEGAZ.*V, 'y'));

WAE = -iterm - jterm-kterm;
WAEa = averageQForcings(WAE, mask, gridvol, ts);

D =  TtoB*GetVar(statefile, diagfile, {'TOTTTEND', 'UDIAG1', '(1)/86400'}, slice);
Dx = Drv(dx, D, 'x');
Dy = Drv(dy, D, 'y');
Dz = Drv(metric, D, 'z');

MC = OMEGAX.*Dx + OMEGAY.*Dy + OMEGAZ.*Dz;
rho = 1031*(1-2e-4.*(THETA - 16.5));
rhodiv = Drv(dx, rho.*U, 'x') + Drv(dy, rho.*V, 'y') + Drv(metric, rho.*W, 'z');

MCa = averageQForcings(MC, mask, gridvol, ts);
MCarho = averageQForcings(Q.*rhodiv./rho, mask, gridvol, ts);
%%
figure('Name', 'Eq 15 balance');
plot(Qa)
hold on
plot(WFa);
plot(WAEa);
plot(MCa);
plot(MCarho);

plot(WFa +WAEa+ MCa+MCarho);
plot(-Fricst - Diast);
% plot(Qda)
hold off

legend('Q', 'WF', 'WAE', 'MC', 'MC_{\rho}', 'Sum', 'Sum J vec');

%%
% Try two ways of calculating Jf

JFz = bx.*Fy - by.*Fx;
JFza = -squeeze(nansum(nansum(FricDiv(:,:,1,:))));
JFzalt = -bx.*dFydx - by.*dFxdy;
JFzalta = -squeeze(nansum(nansum(JFzalt(:,:,1,:))));
scatter(JFza, JFzalta);
% figure
% plot(JFza);
% hold on
% plot(JFzalta);
% hold off
%%
% Surface PV (Schneider notation)
S = OMEGAZ.*b;
Sa = squeeze(nansum(nansum(S(:,:,1,:))));
Sa = Sa.*gridvol./vol;
Sa = Sa(:) - Sa(1);

%%
% Check for importance of cumtrapz vs nansum
FricsS = squeeze(nansum(nansum(JFz(:,yl:end-yl,zl,:).*mask(:,yl:end-yl,zl,:)))).*dx.*dy;
