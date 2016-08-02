function [Q,Qdir, JADVx, JADVy, JADVz, JFx, JFy, JFz, JBx, JBy, JBz, JFzN, JBzN, JFzH, JBzH] = calcQBudgetD(diagfile, statefile, etanfile,extrafile, sizes, slice, dx, dy,dz )
nx = sizes(1); ny = sizes(2); nz=sizes(3); nt=sizes(4);
sliceEta = {slice{1}, slice{2}, [1 1], slice{4}};
TtoB = 9.81.*2e-4;

U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);
% disp(slice{1}); disp(slice{2}); disp(slice{3}); disp(slice{4});
ztmp = ncread(statefile, 'Z');
ztmp = ztmp(1:nz);
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);
zL = length(ztmp);
disp('WARNING: ETAN terms unknown units');
dpdx = GetVar(statefile,diagfile, {'Um_dPHdx','(1)'},slice);
% dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.81*(1)'},sliceEta));
dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN','-(1)/1035'},sliceEta));

dEtadx = DPeriodic(dEtadx, dx,'x');
dpdx = permute(repmat((dEtadx), [1, 1, 1, zL]), [1 2 4 3]) + dpdx;

dpdy = GetVar(statefile,diagfile, {'Vm_dPHdy','(1)'},slice);
% dEtady = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.81*(1)'},sliceEta));  
dEtady = squeeze(GetVar(statefile, etanfile, {'ETAN','-(1)/1035'},sliceEta));  

dEtady = DPeriodic(dEtady, dy, 'y');
% dEtady(:,end-1,:) = 0.5.*dEtady(:,end-2,:);
dpdy = permute(repmat((dEtady), [1, 1, 1, zL]), [1 2 4 3]) + dpdy;

%%
% Calculate buoyancy gradients
% disp('Calculate B gradients')
% bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
% bx(end,:,:,:) = 0.5.*(bx(end-1,:,:,:));
% bx(1,:,:,:) = 0.5.*(bx(2,:,:,:));
b = GetVar(statefile, diagfile, {'b', '(1)'}, slice);
bx = DPeriodic(b, dx,'x');

% by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
by = DPeriodic(b, dy, 'y');
% by(:,1,:,:) = 0;
% by(:,end-1:end,:,:) = 0;
% by(:,end-1,:,:) = 0.5*by(:,end-2,:,:);
% by(:,1,:,:) = 0.5.*by(:,2,:,:);
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);
% bz(:,:,end,:)= 0;
%% 
% Calculate Absolute vorticity terms
% disp('Calculate Abs Vorticity Terms');

%ignore w derivatives for consistency with hydrostatic approx.
%OMEGAX = GetVar(statefile, diagfile, {'WVEL', 'VVEL', 'Dy(1) - Dz(2)'}, slice);
OMEGAX = GetVar(statefile, diagfile, {'VVEL', ' - Dz(1)'}, slice);

%OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'WVEL', 'Dz(1) - Dx(2)'}, slice);
OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);

% OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL', '(1) + Dx(2) - Dy(3)'}, slice);
% V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
OMEGAZ = DPeriodic(V, dx,'x');
% OMEGAZ = GetVar(statefile, diagfile, {'VVEL', 'Dx(1)'},slice);
% Uy = Drv(dy, U, 'y');
% Uy(:,end-1,:,:) = 0.5.*Uy(:,end-2,:,:);
Uy = DPeriodic(U, dy, 'y');
OMEGAZ = OMEGAZ  -Uy;
OMEGAZ = OMEGAZ + 1e-4;
Q = OMEGAX.*bx + OMEGAY.*by + OMEGAZ.*bz; %A more direct definition of Q.
% Q = OMEGAZ.*bz;
%%
%Calculate Q from the momentum budget:
advterms = true;
directADV =true;
diagOUT = false;
LHSUt = GetVar(statefile, diagfile, {'TOTUTEND', '(1)/86400'}, slice);
LHSVt = GetVar(statefile, diagfile, {'TOTVTEND', '(1)/86400'}, slice);
if advterms
    disp('Including Adv Terms in QDir');

    if ~directADV
      VADVTERM = GetVar(statefile, diagfile', { 'Vm_Advec', '(1)'}, slice);
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
        Vx = DPeriodic(V, dx, 'x');
%         Vy = GetVar(statefile, diagfile, {'VVEL', 'Dy(1)'}, slice);
        Vy = DPeriodic(V, dy, 'y');
        Vz = GetVar(statefile, diagfile, {'VVEL', 'Dz(1)'}, slice); % This needs to be fixed
        VADVTERM = -U.*Vx - V.*Vy - W.*Vz + Vm_Cori;
        end
    end
    %Notice I've removed the dpdy here and below...doesn't matter, perfect cancellation (numerics aside).
    LHSV = LHSVt - VADVTERM; - dpdy;
    if ~directADV
    UADVTERM = GetVar(statefile, diagfile, {'Um_Advec', '(1)'}, slice);
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
        Ux = DPeriodic(U, dx,'x');
%         Uy = GetVar(statefile, diagfile, {'UVEL', 'Dy(1)'}, slice);
        Uy = DPeriodic(U, dy, 'y');
        Uz = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);
        UADVTERM = -U.*Ux - V.*Uy - W.*Uz + Um_Cori;
        end

    end
    LHSU = LHSUt - UADVTERM; - dpdx;
   
end

%  LHSU = -dpdx;
QXz = Drv(metric, LHSUt, 'z',1,1);
% QXy = Drv(dy, LHSUt, 'y', 1,1);
QXy = DPeriodic(LHSUt, dy, 'y');
QVz = Drv(metric, LHSVt, 'z', 1, 1);
QVx = DPeriodic(LHSVt, dx, 'x');
Qmom = -bx.*QVz + by.*QXz + bz.*(QVx - QXy);

if advterms
    if ~directADV
    LHSb = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND',  '(1)/86400'}, slice);
    BADV = TtoB.*GetVar(statefile, diagfile, {'UDIAG1', '(1)'}, slice);
    LHSb = LHSb - BADV;
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
LHSbt = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND', '(1)/86400'}, slice);

% LHBx = Drv(dx, LHSb, 'x', 1, 1);
LHBx = DPeriodic(LHSbt, dx,'x');
% LHBy = Drv(dy, LHSbt, 'y', 1, 1);
LHBy = DPeriodic(LHSbt, dy, 'y');
LHBz = Drv(metric, LHSbt, 'z', 1, 1);
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
Fx = Fx - dpdx;


Fy = LHSV;
Fy = Fy - dpdy;

JFx = -bz.*Fy; JFy = bz.*Fx; JFz = bx.*Fy - by.*Fx;

%Horizontal Friction Terms
FxH = GetVar(statefile, extrafile, {'VISCx_Um', 'VISCy_Um', ['Dx(1)/',num2str(dy*dz),'+Dy(2)/',num2str(dx*dz)]}, slice);
FyH = GetVar(statefile, extrafile, {'VISCx_Vm', 'VISCy_Vm', ['Dx(1)/',num2str(dy*dz),'+Dy(2)/',num2str(dx*dz)]}, slice);

JFzH = bx.*FyH - by.*FxH;

%Numeric Friction Terms
FxN =  GetVar(statefile, diagfile, { 'Um_Advec', '(1)'}, slice)-UADVTERM;
FyN =  GetVar(statefile, diagfile, { 'Vm_Advec', '(1)'}, slice)-VADVTERM;

JFzN = bx.*FyN - by.*FxN;

%DIABATIC TERMS
D = LHSb;

JBx = -OMEGAX.*D; JBy = -OMEGAY.*D; JBz = -OMEGAZ.*D;

% Horizontal Diabatic Terms
DH = TtoB.*GetVar(statefile, extrafile, {'DFxE_TH', 'DFyE_TH', ['Dx(1)/',num2str(dy*dz),'+Dy(2)/',num2str(dx*dz)]}, slice);
JBzH = -OMEGAZ.*DH;

%Numeric Diabatic Terms
DN = -TtoB.*GetVar(statefile, diagfile, {'UDIAG1', '(1)'}, slice)-BADV;

JBzN = -OMEGAZ.*DN;
%%

end