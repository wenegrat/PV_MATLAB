function [Q,Qdir, JADVx, JADVy, JADVz, JFx, JFy, JFz, JBx, JBy, JBz] = calcQBudget(diagfile, statefile, etanfile, sizes, slice, dx, dy )
nx = sizes(1); ny = sizes(2); nt=sizes(3);
sliceEta = {slice{1}, slice{2}, [1 1], slice{4}};
TtoB = 9.81.*2e-4;

U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);
% disp(slice{1}); disp(slice{2}); disp(slice{3}); disp(slice{4});
ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);
zL = length(ztmp);

dpdx = GetVar(statefile,diagfile, {'Um_dPHdx','(1)'},slice);
dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.81*(1)'},sliceEta));
dEtadx = DxPeriodic(dEtadx, dx);
dpdx = permute(repmat((dEtadx), [1, 1, 1, zL]), [1 2 4 3]) + dpdx;

dpdy = GetVar(statefile,diagfile, {'Vm_dPHdy','(1)'},slice);
dEtady = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.81*Dy(1)'},sliceEta));  
dpdy = permute(repmat((dEtady), [1, 1, 1, zL]), [1 2 4 3]) + dpdy;

%%
% Calculate buoyancy gradients
% disp('Calculate B gradients')
% bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
% bx(end,:,:,:) = 0.5.*(bx(end-1,:,:,:));
% bx(1,:,:,:) = 0.5.*(bx(2,:,:,:));
b = GetVar(statefile, diagfile, {'b', '(1)'}, slice);
bx = DxPeriodic(b, dx);

by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
by(:,1,:,:) = 0;

by(:,end-1:end,:,:) = 0;
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);
bz(:,:,end,:)= 0;
%% 
% Calculate Absolute vorticity terms
% disp('Calculate Abs Vorticity Terms');

%ignore w derivatives for consistency with hydrostatic approx.
%OMEGAX = GetVar(statefile, diagfile, {'WVEL', 'VVEL', 'Dy(1) - Dz(2)'}, slice);
OMEGAX = GetVar(statefile, diagfile, {'VVEL', ' - Dz(1)'}, slice);

%OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'WVEL', 'Dz(1) - Dx(2)'}, slice);
OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);

% OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL', '(1) + Dx(2) - Dy(3)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
OMEGAZ = DxPeriodic(V, dx);
% OMEGAZ = GetVar(statefile, diagfile, {'VVEL', 'Dx(1)'},slice);
OMEGAZ = OMEGAZ + GetVar(statefile, diagfile, {'UVEL', '-Dy(1)'}, slice);
OMEGAZ = OMEGAZ + 1e-4;
Q = OMEGAX.*bx + OMEGAY.*by + OMEGAZ.*bz; %A more direct definition of Q.
% Q = OMEGAZ.*bz;
%%
%Calculate Q from the momentum budget:
advterms = true;
directADV = false;
diagOUT = true;
LHSU = GetVar(statefile, diagfile, {'TOTUTEND', '(1)/86400'}, slice);
LHSV = GetVar(statefile, diagfile, {'TOTVTEND', '(1)/86400'}, slice);
if advterms
    disp('Including Adv Terms in QDir');

    if ~directADV
      VADVTERM = GetVar(statefile, 'extra.nc', { 'Vm_Advec', '(1)'}, slice);
%     VADVTERM = GetVar(statefile, diagfile, {'ADVx_Vm', 'ADVy_Vm', 'ADVrE_Vm', 'Vm_Cori', ['-Dx(1)/',divstrz, '-Dy(2)/', divstrz,'-Dz(3)/', divstrh,'+(4)']}, slice);
%      VADVTERM = GetVar(statefile, diagfile, {'ADVx_Vm', 'ADVy_Vm', 'ADVrE_Vm',  ['-Dx(1)/',divstrz, '-Dy(2)/', divstrz,'-Dz(3)/', divstrh]}, slice);
%     VADVTERM = GetVar(statefile, diagfile, {'Vm_Cori', '(1)'}, slice);
    else 
        Vm_Cori = GetVar(statefile, diagfile, {'Vm_Cori', '(1)'}, slice);
        if diagOUT
            disp('Using Model ADV diagnostics');
        VADVTERM = -GetVar(statefile,diagfile,{'VADVTERM', '(1)'}, slice);
        VADVTERM = VADVTERM + Vm_Cori;
        else
         disp('Using Direct Advection Terms');
%         Vx = GetVar(statefile, diagfile, {'VVEL', 'Dx(1)'}, slice);
        V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
        Vx = DxPeriodic(V, dx);
        Vy = GetVar(statefile, diagfile, {'VVEL', 'Dy(1)'}, slice);
        Vz = GetVar(statefile, diagfile, {'VVEL', 'Dz(1)'}, slice);
        VADVTERM = -U.*Vx - V.*Vy - W.*Vz + Vm_Cori;
        end
    end
    %Notice I've removed the dpdy here and below...doesn't matter, perfect cancellation (numerics aside).
    LHSV = LHSV - VADVTERM; - dpdy;
    if ~directADV
    UADVTERM = GetVar(statefile, 'extra.nc', {'Um_Advec', '(1)'}, slice);
%      UADVTERM = GetVar(statefile, diagfile, {'ADVx_Um', 'ADVy_Um', 'ADVrE_Um', 'Um_Cori', ['-Dx(1)/',divstrz, '-Dy(2)/', divstrz,'-Dz(3)/', divstrh,'+(4)']}, slice);
%      UADVTERM = GetVar(statefile, diagfile, {'ADVx_Um', 'ADVy_Um', 'ADVrE_Um', ['-Dx(1)/',divstrz, '-Dy(2)/', divstrz,'-Dz(3)/', divstrh]}, slice);
 %    UADVTERM = GetVar(statefile, diagfile, {'Um_Cori', '(1)'}, slice);
   else 
        Um_Cori = GetVar(statefile, diagfile, {'Um_Cori', '(1)'}, slice);
        if diagOUT
            UADVTERM = -GetVar(statefile, diagfile, {'UADVTERM', '(1)'}, slice);
            UADVTERM = UADVTERM + Um_Cori; %On RHS
        else
%         Ux = GetVar(statefile, diagfile, {'UVEL', 'Dx(1)'}, slice);
        U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
        Ux = DxPeriodic(U, dx);
        Uy = GetVar(statefile, diagfile, {'UVEL', 'Dy(1)'}, slice);
        Uz = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);
        UADVTERM = -U.*Ux - V.*Uy - W.*Uz + Um_Cori;
        end

    end
    LHSU = LHSU - UADVTERM; - dpdx;
   
end

%  LHSU = -dpdx;
QXz = Drv(metric, LHSU, 'z',1,1);
QXy = Drv(dy, LHSU, 'y', 1,1);
QVz = Drv(metric, LHSV, 'z', 1, 1);
% QVx = Drv(dx, LHSV, 'x',1,1);
QVx = DxPeriodic(LHSV, dx);
Qmom = -bx.*QVz + by.*QXz + bz.*(QVx - QXy);

if advterms
    if ~directADV
    LHSb = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND',  '(1)/86400'}, slice);
    LHSb = LHSb - TtoB.*GetVar(statefile, 'extra.nc', {'UDIAG1', '(1)'}, slice);
    else
    LHSb = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND', '(1)/86400'}, slice);
        if diagOUT
            BADV = GetVar(statefile, diagfile, {'BADVTERM', '(1)'}, slice);
        else
            BADV = U.*bx + V.*by + W.*bz;
        end
        LHSb = LHSb+BADV;   
    end
else
    LHSb = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND', '(1)/86400'}, slice);
end
%disp('Badv=0');
% LHSb = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND', '(1)/86400'}, slice); %NOTE THAT I AM AD HOC JUST TURNING OFF ADVEC

% LHBx = Drv(dx, LHSb, 'x', 1, 1);
LHBx = DxPeriodic(LHSb, dx);
LHBy = Drv(dy, LHSb, 'y', 1, 1);
LHBz = Drv(metric, LHSb, 'z', 1, 1);
Qb = OMEGAX.*LHBx + OMEGAY.*LHBy + OMEGAZ.*LHBz;
Qdir = Qmom + Qb;

% Qdir = bz.*(QVx-QXy) + OMEGAZ.*LHBz; %Vert Component Only


% 
% %ADVECTIVE TERMS
JADVx = U.*Q; 
JADVz = W.*Q;
JADVy = V.*Q;
%FRICTION TERMS
Fx = LHSU;
% if ~directADV
%     Fx = LHSU;
% %     Fx = GetVar(statefile, diagfile, {'TOTUTEND','Um_Advec', ' (1)/86400 - (2)'},slice);
% %  Fx = GetVar(statefile, diagfile, {'TOTUTEND','Um_Advec','Um_Diss', ' (1)/86400 - (2)+(3)'},slice);
% else
%    Fx = GetVar(statefile, diagfile, {'TOTUTEND', ' (1)/86400 '},slice);
%    UADVTERM = U.*Ux + V.*Uy + W.*Uz - Um_Cori;
%    Fx = Fx + UADVTERM ;
% end
Fx = Fx - dpdx;
% if ~directADV
%     Fy = GetVar(statefile, diagfile, {'TOTVTEND', 'Vm_Advec', '(1)/86400 - (2)'}, slice);
% % Fy = GetVar(statefile, diagfile, {'TOTVTEND', 'Vm_Advec','Vm_Diss', '(1)/86400 - (2)+(3)'}, slice);
% else
%    Fy = GetVar(statefile, diagfile, {'TOTVTEND', ' (1)/86400 '},slice);
%    VADVTERM = U.*Vx + V.*Vy + W.*Vz - Vm_Cori;
%    Fy = Fy + VADVTERM ;
% end
Fy = LHSV;
Fy = Fy - dpdy;

JFx = -bz.*Fy; JFy = bz.*Fx; JFz = bx.*Fy - by.*Fx;

%DIABATIC TERMS
D = LHSb;
% D =  TtoB*GetVar(statefile, diagfile, {'TOTTTEND', 'UDIAG1', '(1)/86400 - (2)'}, slice);
%Using built in diags
%  D =  TtoB*GetVar(statefile, diagfile, {'TOTTTEND', 'UDIAG1', '(1)/86400'}, slice);
%Using the built in flux diagnostics
% ADVx_TH = GetVar(statefile, diagfile, {'ADVx_TH', ['Dx(1)/',divstrz]}, slice);
% ADVy_TH = GetVar(statefile, diagfile, {'ADVy_TH', ['Dy(1)/',divstrz]}, slice);
% ADVr_TH = GetVar(statefile, diagfile, {'ADVr_TH', ['Dz(1)/',divstrh]}, slice);
% ADV= -ADVx_TH-ADVy_TH-ADVr_TH;
% D = D - TtoB.*ADV;

%Hardcoding stuff for testing
% DIFF = GetVar(statefile, diagfile, {'DFxE_TH', 'DFyE_TH','DFrI_TH', '-Dx(1)/2.5e3 - Dy(2)/2.5e3-Dz(3)/1e6'}, slice);
% TFLUX = GetVar(statefile, etanfile, {'TFLUX', '(1)'}, {0, 0, [1 1], slice{4}});
% Cw = 3994;		  
% H = 2.5;
% DIFF(:,:,1,:) = DIFF(:,:,1,:) + TFLUX./(1035*Cw*H);
% D = TtoB.*(DIFF);

JBx = -OMEGAX.*D; JBy = -OMEGAY.*D; JBz = -OMEGAZ.*D;
%%

end