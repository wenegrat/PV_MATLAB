%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL VOlUME ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isoT = [0 100]; %Needed for plot below
mask = ones(nx, ny, nz, tslice(end)-tslice(1)+1);
vol = squeeze(sum(sum(sum(mask.*gridvol))));

%INTEGRATE Q
IntegrateQTerms;

%AreaIntegrateJTerms
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JFb, ~    ] = areaIntegrateJVecs(squeeze(JFz(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JFa = JFs-JFb;
[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JBb, ~    ] = areaIntegrateJVecs(squeeze(JBz(:,:,end,:)), squeeze(mask(:,:,end,:)), dx*dy, ts, vol);
JBa = JBs-JBb;

[JBsz, dJBzdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:) - JBzH(:,:,2,:) - JBzN(:,:,2,:)),squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
%PLOT (QBudget, dQdt)
titleString = ['Full Volume           Surface B_0: ', num2str(squeeze(Q0(1,1,1)))];
QBudgetPlot;
dQdtPlot;

if (saveflag)
    FigString = [IDString, '_FullVol'];
    saveas(QBudgFig, ['QBudget_', FigString, '.png']);
    saveas(dQdtFig, ['dQdt_', FigString, '.png']);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISOPYCNAL ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isoT = [15.0 16.25]; [16.5 16.75];
mask = (THETA(:,:,:,:)>isoT(1)) & (THETA(:,:,:,:)<isoT(2));
vol = squeeze(sum(sum(sum(mask.*gridvol))));

%INTEGRATE Q
IntegrateQTerms;

%AreaIntegrateJTerms
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol); %Surface Only by impermeability
JFa = JFs;
[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
JBa = JBs;

%PLOT (QBudget, dQdt)
titleString = ['Iso Volume           Surface B_0: ', num2str(squeeze(Q0(1,1,1)))];
QBudgetPlot;
dQdtPlot;

if (saveflag)
    FigString = [IDString, '_ISO', num2str(isoT(1)),'-', num2str(isoT(2))];
    saveas(QBudgFig, ['QBudget_', FigString, '.png']);
    saveas(dQdtFig, ['dQdt_', FigString, '.png']);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAYER ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isoT = [1 100];
mask = zeros(nx, ny, nz, nt);
depthind = 35;
mask(:,:,1:depthind,:) = 1; %z > fixed depth

% deltaCrit = 0.01./(-2e-4*1035); % Convert to a delta T criteria.
% tprime = THETA - repmat(THETA(:,:,1,:), [1, 1, nz, 1]);
% mask = tprime>deltaCrit;
vol = squeeze(sum(sum(sum(mask.*gridvol))));
% 
%INTEGRATE Q
IntegrateQTerms;
% 
% %Volume Integrals of J divergences
% areaIntegrateJDiv;

%AreaIntegrateJTerms
[JFs, dJFdt] = areaIntegrateJVecs(squeeze(JFz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JFb, ~    ] = areaIntegrateJVecs(squeeze(JFz(:,:,depthind,:)), squeeze(mask(:,:,depthind,:)), dx*dy, ts, vol);
Fric = JFs-JFb;
[JBs, dJBdt] = areaIntegrateJVecs(squeeze(JBz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JBb, ~    ] = areaIntegrateJVecs(squeeze(JBz(:,:,depthind,:)), squeeze(mask(:,:,depthind,:)), dx*dy, ts, vol);
Dia = JBs-JBb;
[JAs, dJAdt] = areaIntegrateJVecs(squeeze(JAz(:,:,2,:)), squeeze(mask(:,:,2,:)), dx*dy, ts, vol);
[JAb, ~    ] = areaIntegrateJVecs(squeeze(JAz(:,:,depthind,:)), squeeze(mask(:,:,depthind,:)), dx*dy, ts, vol);
Adv = JAs-JAb;
% Adv = -JAb;
%PLOT (QBudget, dQdt)
titleString = ['Layer Analysis           Surface B_0: ', num2str(squeeze(Q0(1,1,1)))];
QBudgetPlotAdv;
% dQdtPlot;
if (saveflag)
    FigString = [IDString, '_Layer_', num2str(isoT(2))];
    saveas(QBudgAdvFig, ['QBudget_', FigString, '.png']);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confirm Numerics and Horizontal Terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConfirmNumHoriz;
if (saveflag)
    FigString = [IDString];
    saveas(NumHorizFig, ['NumHorizTerms_', FigString, '.png']);
end

%%
RiPlot;
if (saveflag)
    FigString = [IDString];
    saveas(RiPlotFig, ['RiPlot_', FigString, '.png']);
end
%%
toc./60 % in minutes