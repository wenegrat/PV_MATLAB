%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Balance ZONAL Mom Budget as test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('TOTUTEND');
statefile = 'state.nc';
diagfile = 'diag.nc';
etanfile = 'etan.nc';
TOTUTEND = GetVar(statefile,diagfile, {'TOTUTEND', '(1)/86400'},{0, 0, 0, 0});
[nx ny nz nt] = size(TOTUTEND);
disp('RHS');
% RHS = GetVar('state_0000000000.nc','fullDiag_0000000000.nc', {'Um_dPHdx','Um_Diss', 'Um_Ext', 'VISrI_Um'...
%       '(1) + (2) + (3) +Dz(4)./1e6'},{0, 0, 0, 0});
dpdx = GetVar(statefile,diagfile, {'Um_dPHdx','(1)'},{0, 0, 0, 0});
Eta = GetVar(statefile, etanfile, {'ETAN', '(1)'}, {0, 0, [1 1], 0});
dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN',...
      '-9.8*Dx(1)'},{0, 0,[1 1], 0}));
dpdx = permute(repmat(squeeze(dEtadx), [1, 1, 1, nz]), [1 2 4 3]) + dpdx;

diss = GetVar(statefile,diagfile, {'Um_Diss','(1)'},{0, 0, 0, 0});
advec = GetVar(statefile,diagfile, {'Um_Advec','(1)'},{0, 0, 0, 0});
visi = GetVar(statefile,diagfile, {'VISrI_Um','- Dz(1)/1e6'},{0, 0, 0, 0});
% abgu = GetVar(statefile,diagfile, {'AB_gU','(1)'},{0, 0, 0, 0});

RHSufull = dpdx + diss +advec+visi;%+abgu;
% RHSu = GetVar('state_0000000000.nc','fullDiag_0000000000.nc', {'Um_dPHdx','Um_Diss', 'Um_Advec', 'VISrI_Um', 'AB_gU'...
%       '(1) + (2) + (3) - Dz(4)/1e6 + (5)'},{0, 0, 0, 0});
% disp('dEta');
% 
% dEtadx = squeeze(GetVar('state_0000000000.nc','etanDiag_0000000000.nc', {'ETAN',...
%       '-9.8*Dx(1)'},{0, 0,[1 1], 0}));

% 
% RHSufull = permute(repmat(squeeze(dEtadx),[1, 1, 1, nz]), [1 2 4 3]) + RHSu;

%%
indj = 20; indk =20; d = 1;

plot(squeeze(TOTUTEND(indj, indk, d,:)), 'LineWidth', 2);
hold on
% plot(squeeze(dpdx(indj, indk, d,:)));
% plot(squeeze(diss(indj, indk, d,:)));
% plot(squeeze(advec(indj, indk, d,:)));
% plot(squeeze(visi(indj, indk, d,:)));
% plot(squeeze(abgu(indj, indk, d,:)));
plot(squeeze(RHSufull(indj, indk, d,:)));
hold off
legend('dU/dt', 'dPdx', 'Diss', 'Advec', 'Vis_I', 'AB', 'RHS');
title('Zonal Mom Budget')
ylabel('m/s^2'); xlabel('Hours');
grid on
%%
% Plot terms.
indj = 10; indk =10; d = 80;

subplot(2,1,1)
scatter(squeeze(TOTUTEND(indj, indk, d,:)), squeeze(RHSufull(indj, indk, d, :)));
grid on
subplot(2,1,2)
plot(squeeze(TOTUTEND(indj, indk, d,:)));
hold on
plot( squeeze(RHSufull(indj, indk, d, :)));
hold off
grid on
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Balance MERIDIONAL Mom Budget as test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('TOTVTEND');
TOTVTEND = GetVar('state_0000000000.nc', 'fullDiag_0000000000.nc', {'TOTVTEND', '(1)/86400'},{0, 0, 0, 0});
disp('RHS');
% RHS = GetVar('state_0000000000.nc','fullDiag_0000000000.nc', {'Um_dPHdx','Um_Diss', 'Um_Ext', 'VISrI_Um'...
%       '(1) + (2) + (3) +Dz(4)./1e6'},{0, 0, 0, 0});
RHSv = GetVar('state_0000000000.nc','fullDiag_0000000000.nc', {'Vm_dPHdy','Vm_Diss', 'Vm_Advec', 'VISrI_Vm', 'AB_gV'...
      '(1) + (2) + (3) - Dz(4)/1e6 + (5)'},{0, 0, 0, 0});
disp('dEta');

dEtady = squeeze(GetVar('state_0000000000.nc','etanDiag_0000000000.nc', {'ETAN',...
      '-9.8*Dy(1)'},{0, 0,[1 1], 0}));
RHSvfull = permute(repmat(squeeze(dEtady),[1, 1, 1, nz]), [1 2 4 3]) + RHSv;
%%
% Plot terms.
indj = 30; indk =10; d = 10;

subplot(2,1,1)
scatter(squeeze(TOTVTEND(indj, indk, d,:)), squeeze(RHSvfull(indj, indk, d, :)));
grid on
subplot(2,1,2)
plot(squeeze(TOTVTEND(indj, indk, d,:)));
hold on
plot( squeeze(RHSvfull(indj, indk, d, :)));
hold off
grid on
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Balance Temp Budget as test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
statefile = 'state.nc';'state_0000000000.nc';
diagfile =  'diag.nc'; 'fullDiag_0000000000.nc';
etanfile = 'etan.nc';
TOTTTEND = GetVar(statefile,diagfile, {'TOTTTEND', '(1)/86400'},{0, 0, 0, 0});
% ADVT = GetVar('state_0000000000.nc', 'fullDiag_0000000000.nc', {'ADVr_TH', 'ADVx_TH', 'ADVy_TH',...
%      '-Dz(1)/1e6 + Dx(2)/5e3 + Dy(3)/5e3'},{0, 0, 0, 0});
% DIFT = GetVar('state_0000000000.nc', 'fullDiag_0000000000.nc', {'DFrE_TH', 'DFxE_TH', 'DFyE_TH','DFrI_TH',...
%      '-Dz(1)/1e6 + Dx(2)/5e3 + Dy(3)/5e3 - Dz(4)/1e6'},{0, 0, 0, 0});
KPPT = GetVar(statefile, diagfile, {'KPPg_TH',...
     '-Dz(1)/1e6'},{0, 0, 0, 0});
THETA = GetVar(statefile, diagfile, {'THETA', '(1)'},{0, 0, 0, 0});

ADVD = GetVar(statefile, diagfile, {'UDIAG1', '(1)'}, {0, 0, 0, 0});
% ADVTZ = GetVar(statefile, diagfile, {'ADVr_TH', '-Dz(1)/1e6'},{0, 0, 0, 0});
% ADVTX = GetVar(statefile, diagfile, {'ADVx_TH', 'Dx(1)/5e3'},{0, 0, 0, 0});
% ADVTY = GetVar(statefile, diagfile, {'ADVy_TH', 'Dy(1)/5e3'},{0, 0, 0, 0});

DIFTZ = GetVar(statefile, diagfile, {'DFrE_TH', '-Dz(1)/1e6'},{0, 0, 0, 0});
DIFTIZ = GetVar(statefile, diagfile, {'DFrI_TH', '-Dz(1)./1e6'},{0, 0, 0, 0});
DFRI_TH = GetVar(statefile, diagfile, {'DFrI_TH', '(1)'}, {0,0,0,0});
DIFTX = GetVar(statefile, diagfile, {'DFxE_TH', 'Dx(1)/5e3'},{0, 0, 0, 0});
DIFTY = GetVar(statefile, diagfile, {'DFyE_TH', 'Dy(1)/5e3'},{0, 0, 0, 0});
% ABGT = GetVar(statefile, diagfile, {'AB_gT', '(1)'},{0, 0, 0, 0});
WTHMASS = GetVar(statefile, diagfile,{'WTHMASS', '(1)'},{0, 0, 0, 0});
%Test nondivergence of velocity derivative operator.
% NONDIV = GetVar(statefile, diagfile, {'UVEL', 'VVEL', 'WVEL', 'Dx(1)+Dy(2) + Dz(3)'}, {0 ,0, 0, 0});
TFLUX = GetVar(statefile, etanfile, {'TFLUX', '(1)'}, {0, 0, [1 1], 0});
[nx ny nd nt] = size(THETA);
TFLUXF = zeros(nx, ny, nd, nt);
TFLUXF(:,:,1,:) = TFLUX(:,:,:, 1:nt);
alpha = 2e-4;
Cw    = 3994;		  
H = 5; % This is a hack, basically putting everything into upper two bins...
TFLUXF = TFLUXF./(1031*Cw*H);
%%
surfCorrTerm = repmat((repmat(nansum(nansum(WTHMASS(:,:,1,:)))./(48*48),[48 48 1 1])  - WTHMASS(:,:,1,:))./(5), [1 1 100 1]);

%
dTdt = NaN(nx, ny, nz, nt);
for i=1:nx;
    for j=1:ny;
        for d = 1:100
            dTdt(i,j,d,:) = gradient(squeeze(THETA(i,j,d,:)), 3600);
        end
    end
end
%%

% Compare ADVD with NONDIVERGENT ADVECTIVE FLUXES

indj = 4; indk = 3, d = 5;
NONDIVA = ADVTZ - ADVTX - ADVTY + NONDIV.*THETA;

scatter(squeeze(ADVD(indj, indk, d, :)), squeeze(NONDIVA(indj, indk, d, :)));
grid on

%%
% Volume integrated budget

LHSA = squeeze(nanmean(nanmean(nanmean(TOTTTEND))));
% LHSA2 = squeeze(nanmean(nanmean(nanmean(dTdt))));

ADVA = squeeze(nanmean(nanmean(nanmean(ADVD))));
TFLUXA = squeeze(nanmean(nanmean(nanmean(TFLUXF))));
DIFFA = squeeze(nanmean(nanmean(nanmean(DIFT))));
KPPA = squeeze(nanmean(nanmean(nanmean(KPPT))));

plot(LHSA);
hold on
plot(ADVA);
plot(TFLUXA);
plot(DIFFA);
plot(KPPA+ADVA+TFLUXA-DIFFA);
plot(KPPA)
% plot(LHSA2, '--');
hold off
legend('LHS', 'ADV', 'TFLUX', 'DIFF', 'SUM', 'KPP');
%%
% Plot terms.
indj = 40; indk =40; d = 1; 
st = 1;
% ADVTZTemp = ADVTZ;
% ADVTZTemp(:,:,1:99,:) = ADVTZTemp(:,:,2:end,:);
% ADVT =  -ADVTZTemp+ADVTX + ADVTY;
DIFT = -DIFTZ - DIFTIZ + DIFTX + DIFTY;

RHST = +KPPT - DIFT + ADVD+ TFLUXF;%- surfCorrTerm/30;%- NONDIV.*THETA;

% subplot(2,1,1)
% scatter(squeeze(TOTTTEND(indj, indk, d,:)), squeeze(RHST(indj, indk, d, :)));
% grid on
% subplot(2,1,2)
plot((squeeze(TOTTTEND(indj, indk, d,:))), 'LineWidth', 2);
hold on
plot(( squeeze(RHST(indj, indk, d, :))), '--', 'LineWidth', 2);
% plot(nanmean(squeeze(RHST(indj, indk, d:d+1, :))), '--', 'LineWidth', 2);

% plot( squeeze(surfCorrTerm(indj, indk, d+1, :)));
% plot( squeeze(ABGT(indj, indk, d, :)));

% plot( squeeze(-ADVTZ(indj, indk, d, :)));
% plot( squeeze(ADVTX(indj, indk, d, :)));
% plot( squeeze(ADVTY(indj, indk, d, :)));

% plot( squeeze(-DIFT(indj, indk, d, :)));
% plot( squeeze(-KPPT(indj, indk, d, :)));
hold off
grid on
legend('TEND', 'RHS', 'ADVT', 'DIFT', 'KPPT');
title('Temperature Budget');
ylabel('^{\circ}/s'); xlabel('Hours');
