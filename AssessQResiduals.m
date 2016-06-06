%%
% Try to assess residuals in the Q budget.
 close all; clc;
 statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
% Parameters            
dx = 1000; dy = dx; dz = 5;
divis = '1e6';
ts = 3600;
TtoB = 9.81.*2e-4;
indj = 1; indk = 20; d = 20;
tslice = [400 480];
slice={0, 0, 0, tslice};%100 120
sliceEta={0,0,[1 1],tslice};%251 271


%ZONAL MOM BUDGET
disp('Calculating Zonal Mom Budget Terms');
TOTUTEND = GetVar(statefile,diagfile, {'TOTUTEND', '(1)/86400'},slice);
[xL, yL, zL, tL] = size(TOTUTEND);
dpdx = GetVar(statefile,diagfile, {'Um_dPHdx','(1)'},slice);
dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.8*Dx(1)'},sliceEta));
if ndims(dEtadx)==2
    disp('NDIMS == 2');
    temp = dEtadx;
    clear dEtadx;
    dEtadx(1,:,:) = zeros(size(temp));
end
dpdx = permute(repmat((dEtadx), [1, 1, 1, zL]), [1 2 4 3]) + dpdx;
% RHSDISSU =  GetVar(statefile,diagfile, {'Um_Diss', 'VISrI_Um','VISCx_Um', 'VISCy_Um','(1)+Dz(2)/2.5e5+Dx(3)/1250 +Dy(4)/1250'},slice);
RHSDISSU =  GetVar(statefile,diagfile, {'Um_Diss', 'VISrI_Um',['(1)+Dz(2)/',divis]},slice);

RHSADVECU = GetVar(statefile, diagfile, {'Um_Advec', '(1)'}, slice);
RHSU = RHSDISSU + RHSADVECU + dpdx;
%Meridional MOM BUDGET
disp('Calculating Meridional Mom Budget Terms');
TOTVTEND = GetVar(statefile,diagfile, {'TOTVTEND', '(1)/86400'},slice);
dpdy = GetVar(statefile,diagfile, {'Vm_dPHdy','(1)'},slice);
dEtady = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.8*Dy(1)'},sliceEta));  
if ndims(dEtady)==2
    disp('NDIMS == 2');
    temp = dEtady;
    clear dEtady;
    dEtady(1,:,:) = temp;
end
dpdy = permute(repmat((dEtady), [1, 1, 1, zL]), [1 2 4 3]) + dpdy;
RHSDISSV =  GetVar(statefile,diagfile, {'Vm_Diss', 'VISrI_Vm',['(1)+Dz(2)/',divis]},slice);
RHSADVECV = GetVar(statefile, diagfile, {'Vm_Advec', '(1)'}, slice);
RHSV = RHSDISSV + RHSADVECV + dpdy;


ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, xL, yL, tL]), [2 3 1 4]);
%%
d=1;
figure 
plot(squeeze(TOTUTEND(indj, indk, d,:)), 'LineWidth', 2);
hold on
plot(squeeze(RHSU(indj, indk, d,:)), 'LineWidth', 2, 'LineStyle', '--');
plot(squeeze(TOTVTEND(indj, indk, d,:)), 'LineWidth', 2);
plot(squeeze(RHSV(indj, indk, d,:)), 'LineWidth', 2, 'LineStyle', '--');
hold off
legend('dU/dt', 'RHS_U', 'dV/dt', 'RHS_V');
%%
%Calc RESIDUALS
disp('Calculate Mom Residuals')
RESU = TOTUTEND - RHSU;
RESV = TOTVTEND - RHSV;
RESUz = Drv(metric, RESU, 'z', 1, 1);
RESUy = Drv(dy, RESU, 'y',1, 1);
RESVz = Drv(metric, RESV, 'z',1, 1);
RESVx = Drv(dx, RESV, 'x',1, 1);
%%
% Calculate buoyancy gradients
disp('Calculate B gradients')
bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
if ndims(squeeze(bx))==3; bx = zeros(size(bx)); end
by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);

%%
% This is the 'residual pv' from the residuals of the momentum equation.
disp('Define Q Mom residuals');
RESQM = -bx.*RESVz-by.*RESUz+bz.*(RESVx-RESUy);

%%
disp('Define Q direct from Mom');
LHSU = TOTUTEND-RHSADVECU;
LHSV = TOTVTEND - RHSADVECV;
QXz = Drv(metric, LHSU, 'z',1,1);
QXy = Drv(dy, LHSU, 'y', 1,1);
QVz = Drv(metric, LHSV, 'z', 1, 1);
QVx = Drv(dx, LHSV, 'x',1,1);
QdirMom = -bx.*QVz - by.*QXz + bz.*(QVx - QXy);
if ndims(squeeze(QVx))==3;
    QdirMom = -bx.*QVz - by.*QXz + bz.*( - QXy);
end
%From taking the curl of the momentum equation, ignoring w terms.

%%
% Buoyancy equation residuals
disp('Calculate B Residuals');
TOTTTEND = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND', '(1)/86400'},slice);
RHST = GetVar(statefile, diagfile, {'KPPg_TH','UDIAG1','DFrI_TH', '-Dz(1)/1e6+(2)-Dz(3)/1e6'}, slice);
RHSTADV = TtoB.*GetVar(statefile, diagfile, {'UDIAG1', '(1)'}, slice);

TFLUX = GetVar(statefile, etanfile, {'TFLUX', '(1)'}, {0, 0, [1 1], slice{4}});
[nx ny nd nt] = size(TOTTTEND);
TFLUXF = zeros(nx, ny, nd, nt);
TFLUXF(:,:,1,:) = TFLUX(:,:,:, 1:nt);
Cw    = 3994;		  
H = 5;
TFLUXF = TFLUXF./(1031*Cw*H);
RHST = RHST + TFLUXF;
RHST = TtoB.*RHST;
%%
figure
plot(squeeze(TOTTTEND(indj, indk, d,:)), 'LineWidth', 2);
hold on
plot(squeeze(RHST(indj, indk, d, :)), 'LineWidth', 2, 'LineStyle', '--');
hold off
title('Buoyancy Budget');

%%
disp('Calculate B residual gradients');
REST = TOTTTEND - RHST;
RESTx = Drv(dx, REST, 'x', 1, 1);
RESTy = Drv(dy, REST, 'y', 1, 1);
RESTz = Drv(metric, REST, 'z', 1, 1);
%% 
% Calculate Absolute vorticity terms
disp('Calculate Abs Vorticity Terms');

%ignore w derivatives for consistency with hydrostatic approx.
OMEGAX = GetVar(statefile, diagfile, {'WVEL', 'VVEL', 'Dy(1) - Dz(2)'}, slice);
OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'WVEL', 'Dz(1) - Dx(2)'}, slice);

% OMEGAX = GetVar(statefile, diagfile, { 'VVEL', ' - Dz(1)'}, slice);
% OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'Dz(1) '}, slice);
OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL', '(1) + Dx(2) - Dy(3)'}, slice);
if ndims(squeeze(OMEGAZ))==3; 
    OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'UVEL', '(1)  - Dy(2)'}, slice);
    OMEGAY = Drv(metric, GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice),'z');
end

RESQB = OMEGAX.*RESTx + OMEGAY.*RESTy+OMEGAZ.*RESTz;
%This is the residual in the PV equation from the buoyancy equation.

%%
display('Calculate B contributions to Qdir')
LHSb = TOTTTEND - RHSTADV;
LHBx = Drv(dx, LHSb, 'x', 1, 1);
LHBy = Drv(dy, LHSb, 'y', 1, 1);
LHBz = Drv(metric, LHSb, 'z', 1, 1);

Qb = OMEGAX.*LHBx + OMEGAY.*LHBy + OMEGAZ.*LHBz;
if ndims(squeeze(OMEGAX))==3
Qb =  OMEGAY.*LHBy + OMEGAZ.*LHBz;
end
Qdir = QdirMom - Qb;
%Note that this is the time derivative of Q, as inferred through the mom
%equation.

%%
TOTALQ = RESQM - RESQB;
indj =1; indk = 20; d=3;
figure
subplot(2,1,1)
plot(squeeze(TOTALQ(indj, indk, d, :)), 'LineWidth', 2);
title('Pointwise remainder of PV Budget');

subplot(2,1,2)
plot(squeeze(nanmean(nanmean(nanmean(TOTALQ)))), 'LineWidth', 2);
title('Volume integral of PV budget remainder');
%%
% Calculate Q
disp('Calculate Q and dQ/dt');
Q = OMEGAX.*bx + OMEGAY.*by + OMEGAZ.*bz; %A more direct definition of Q.
QT = NaN(size(Q));
for i=1:nx
    for j=1:ny
        for k =1:zL
        QT(i,j,k,:) = gradient(squeeze(Q(i,j,k,:)), ts);
        end
    end
end
%%
figure
subplot(2,1,1)
plot((squeeze(nanmean(nanmean(nanmean(Qdir))))), 'LineWidth', 2);
hold on
plot((squeeze(nanmean(nanmean(nanmean(QT))))), 'LineWidth', 2, 'LineStyle', '--');
hold off
legend('dQ/dt from Mom', 'dQ/dt from Omega');
subplot(2,1,2)
plot(cumtrapz(squeeze(nanmean(nanmean(nanmean(Qdir))))*ts), 'LineWidth', 2);
hold on
QInt = (squeeze(nanmean(nanmean(nanmean(Q)))));
plot(QInt-QInt(1), 'LineWidth', 2, 'LineStyle', '--');
hold off
legend('Q from Mom', 'Q from Omega');

%%
% Calculate J Vectors
disp('Calculate J Vectors');
U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

%ADVECTIVE TERMS
JADVx = U.*Q; JADVy = V.*Q; JADVz = W.*Q;

%FRICTION TERMS
Fx = GetVar(statefile, diagfile, {'TOTUTEND','Um_Advec', ' (1)/86400 - (2)'},slice);
Fx = Fx - dpdx;
Fy = GetVar(statefile, diagfile, {'TOTVTEND', 'Vm_Advec', '(1)/86400 - (2)'}, slice);
Fy = Fy - dpdy;

JFx = -bz.*Fy; JFy = bz.*Fx; JFz = bx.*Fy - by.*Fx;

%DIABATIC TERMS
D = TtoB*GetVar(statefile, diagfile, {'TOTTTEND', 'UDIAG1', '(1)/86400 - (2)'}, slice);
JBx = -OMEGAX.*D; JBy = -OMEGAY.*D; JBz = -OMEGAZ.*D;

GRADJ = Drv(dx, JADVx+JFx+JBx, 'x') + Drv(dy, JADVy+JFy+JBy, 'y') + Drv(metric, JADVz+JFz+JBz, 'z');


QRES = QT + GRADJ;

%%
d=3;
figure
subplot(2,1,1)
plot(squeeze(JFx(indj, indk, d, :)), 'LineWidth', 2);
hold on
plot(squeeze(JFy(indj, indk, d, :)), 'LineWidth', 2);
plot(squeeze(JFz(indj, indk, d, :)), 'LineWidth', 2);
hold off
title('Jf');
legend('Jfx','Jfy', 'Jfz');

subplot(2,1,2)
plot(squeeze(JBx(indj, indk, d, :)), 'LineWidth', 2);
hold on
plot(squeeze(JBy(indj, indk, d, :)), 'LineWidth', 2);
plot(squeeze(JBz(indj, indk, d, :)), 'LineWidth', 2);
hold off
title('Jbz')
legend('JBx', 'JBy', 'JBz');
%%
% Checking to make sure the residuals from the J vectors aren't
% substantially larger than the residuals from the mom equations
% themselves.
figure
indj = 1; indk = 20; d = 20;
subplot(2,1,1)
plot(squeeze(QRES(indj, indk, d, :)), 'LineWidth', 2);
hold on
plot(squeeze(TOTALQ(indj, indk, d, :)), 'LineWidth', 2, 'LineStyle', '--');
hold off
legend('QRES', 'TOTALQ');
title('Arbitrary Pointwise Budget');

subplot(2,1,2)
plot(squeeze(nanmean(nanmean(nanmean(Qdir)))), 'LineWidth', 2);
hold on
plot(squeeze(nanmean(nanmean(nanmean(GRADJ)))), 'LineWidth', 2, 'LineStyle', '--');
hold off
legend('QRES', 'TOTALQ');
title('Volume Averaged');

%%
disp('Try Volume Integral');
nx = 48; ny = 48;

area = nx.*dx.*ny.*dy;
%Average Horizontally
if ndims(squeeze(JFz))==4
Jbah = squeeze(nansum(nansum(JBz(:,:,:,:))))*dx*dy./area;
Jfah = squeeze(nansum(nansum(JFz(:,:,:,:))))*dx*dy./area;
qah = squeeze(nansum(nansum(Q))).*dx.*dy./area;
qdirah = squeeze(nansum(nansum(Qdir))).*dx.*dy./area;

else %2 Dimensional
    Jbah = squeeze(nansum(JBz)).*dy./(ny.*dy);
    Jfah = squeeze(nansum(JFz)).*dy./(ny.*dy);
    qah = squeeze(nansum(Q)).*dy./(ny.*dy);
    qdirah  = squeeze(nansum(Qdir)).*dy./(ny*dy);
end
Jfahy = squeeze(nansum(JFy(:, end-1, :,:)-JFy(:,2,:,:))).*dx./area;
Jfay = nansum(Jfahy).*dz;
Jfay = cumtrapz(Jfay).*ts;

%Integrating the vertical difference in time.
Jbza = cumtrapz(Jbah(1,:) - Jbah(end,:)).*ts;
% Jbza = cumtrapz(Jbah(2,:)).*ts;

Jfza = cumtrapz(Jfah(1,:) - Jfah(end,:)).*ts;

% Jfza = cumtrapz(Jfah(2,:)).*ts;

%Integrate vertically
qa = nansum(qah)*dz;
qa = qa - qa(1); %Define delta quantity
% qa = -nansum(qah)*dz; %integrate in z
% qa = cumtrapz(qa)*ts; %integrate in time.

qdira = -nansum(qdirah)*dz; %integrate in z
qdira = cumtrapz(qdira)*ts; %integrate in time.
%%
figure
plot(qa, 'LineWidth', 2)
hold on
plot(-Jbza, 'LineWidth', 2);
plot(-Jfza, 'LineWidth', 2);

plot(-(Jbza+Jfza), 'LineWidth', 2, 'LineStyle', '--');
% plot(qdira);
plot(-Jfay);
hold off
legend('Q', 'J_B^z', 'J_F^z', '-SUM');
xlabel('Num time steps (Hourly)');
ylabel('\Delta Q');
grid on

% %%
% scatter(qa + (Jbza + Jfza), 100*Jfay)
% grid on