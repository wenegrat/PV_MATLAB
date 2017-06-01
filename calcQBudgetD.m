function [Q, JADVx, JADVy, JADVz, JFx, JFy, JFz, JBx, JBy, JBz, JFzN, JBzN, JFzH, JBzH, OMEGAZs] = calcQBudgetD(diagfile, statefile, etanfile,extrafile, sizes, slice, dx, dy,dz )
nx = sizes(1); ny = sizes(2); nz=sizes(3); nt=sizes(4);
sliceEta = {slice{1}, slice{2}, [1 1], slice{4}};
TtoB = 9.81.*2e-4;

% Load Velocities
U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

ztmp = ncread(statefile, 'Z');
ztmp = ztmp(1:nz);
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);
zL = length(ztmp);

% Calculate Pressure Gradients
dpdx = GetVar(statefile,diagfile, {'Um_dPHdx','(1)'},slice);
dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.81*(1)'},sliceEta));
dEtadx = DPeriodic(dEtadx, dx,'x');
dpdx = permute(repmat((dEtadx), [1, 1, 1, zL]), [1 2 4 3]) + dpdx;

dpdy = GetVar(statefile,diagfile, {'Vm_dPHdy','(1)'},slice);
dEtady = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.81*(1)'},sliceEta));  
dEtady = DPeriodic(dEtady, dy, 'y');
dpdy = permute(repmat((dEtady), [1, 1, 1, zL]), [1 2 4 3]) + dpdy;

% Calculate buoyancy gradients
b = GetVar(statefile, diagfile, {'b', '(1)'}, slice);
bx = DPeriodic(b, dx,'x');

by = DPeriodic(b, dy, 'y');
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);
% bz(:,:,end,:)= 0;
%% 
% Calculate Absolute vorticity terms

%ignore w derivatives for consistency with hydrostatic approx.
%OMEGAX = GetVar(statefile, diagfile, {'WVEL', 'VVEL', 'Dy(1) - Dz(2)'}, slice);
Wy = DPeriodic(W, dy, 'y');
Wx = DPeriodic(W, dx, 'x');


OMEGAX = GetVar(statefile, diagfile, {'VVEL', ' - Dz(1)'}, slice);
OMEGAX = Wy + OMEGAX;
%OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'WVEL', 'Dz(1) - Dx(2)'}, slice);
OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);
OMEGAY = -Wx + OMEGAY;

% OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL', '(1) + Dx(2) - Dy(3)'}, slice);
% V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
OMEGAZ = DPeriodic(V, dx,'x'); 
Uy = DPeriodic(U, dy, 'y');
OMEGAZ = OMEGAZ  - Uy;
OMEGAZ = OMEGAZ + 1e-4;
OMEGAZs = squeeze(OMEGAZ(:,:,2));
Q = OMEGAX.*bx + OMEGAY.*by + OMEGAZ.*bz; %A more direct definition of Q.

%%
%Calculate Nonconservative Terms from Mom Equation
%%%%%%%%%%%%%%%%%%%%%%%
% Best way is advterms = true; directADV = true; diagOUT = false;

advterms = true; % This flag turns on the advection terms in the momentum equation for flux terms
directADV = true; % This flag is for calculating offline advective terms (isolating numerical non-conservative)
diagOUT = false; % If True use model diagnostics, and add Cori term (depends on mom advection scheme).

LHSUt = GetVar(statefile, diagfile, {'TOTUTEND', '(1)/86400'}, slice);
LHSVt = GetVar(statefile, diagfile, {'TOTVTEND', '(1)/86400'}, slice);

if advterms
    disp('Including Adv Terms in QDir');

    if ~directADV
      VADVTERM = GetVar(statefile, diagfile', { 'Vm_Advec', '(1)'}, slice);
    else 
        Vm_Cori = GetVar(statefile, diagfile, {'Vm_Cori', '(1)'}, slice);
        if diagOUT
            disp('Using Model ADV diagnostics');
        VADVTERM = -GetVar(statefile,diagfile,{'VADVTERM', '(1)'}, slice);
        VADVTERM = VADVTERM + Vm_Cori;
        else
         disp('Using Direct Advection Terms');
        Vx = DPeriodic(V, dx, 'x');
        Vy = DPeriodic(V, dy, 'y');
        Vz = Drv(metric, V, 'z');
        VADVTERM = -U.*Vx - V.*Vy - W.*Vz + Vm_Cori; % Sign convenction is as a RHS term
        end
    end
    %Notice I've removed the dpdy here and below...doesn't matter, perfect cancellation (numerics aside).
    LHSV = LHSVt - VADVTERM;% - dpdy;
    if ~directADV
    UADVTERM = GetVar(statefile, diagfile, {'Um_Advec', '(1)'}, slice);
   else 
        Um_Cori = GetVar(statefile, diagfile, {'Um_Cori', '(1)'}, slice);
        if diagOUT
            UADVTERM = -GetVar(statefile, diagfile, {'UADVTERM', '(1)'}, slice);
            UADVTERM = UADVTERM + Um_Cori; %On RHS
        else
        Ux = DPeriodic(U, dx,'x');
        Uy = DPeriodic(U, dy, 'y');
        Uz = Drv(metric, U, 'z');
        UADVTERM = -U.*Ux - V.*Uy - W.*Uz + Um_Cori;
        end

    end
    LHSU = LHSUt - UADVTERM;% - dpdx;
   
end

% Buoyancy Non-Conservative Terms
if advterms
    if ~directADV
    LHSb = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND',  '(1)/86400'}, slice);
    BADV = -TtoB.*GetVar(statefile, diagfile, {'UDIAG1', '(1)'}, slice);
    LHSb = LHSb + BADV;
    else
    LHSbt = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND', '(1)/86400'}, slice);
        if diagOUT
            BADV = GetVar(statefile, diagfile, {'BADVTERM', '(1)'}, slice);
        else
            BADV = U.*bx + V.*by + W.*bz;
        end
        LHSb = LHSbt+BADV;   % Defined as a LHS term
    end
else
    LHSb = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND', '(1)/86400'}, slice);
end

%%% CALCULATE ALL J Vectors

% ADVECTIVE J Vectors
JADVx = U.*Q; 
JADVz = W.*Q;
JADVy = V.*Q;

% FRICTION J Vectors
Fx = LHSU;
Fx = Fx - dpdx;
Fy = LHSV;
Fy = Fy - dpdy;

JFx = -bz.*Fy; JFy = bz.*Fx; JFz = bx.*Fy - by.*Fx;

% Isolate Horizontal Friction Terms
% Note this is not for the budget, but just for sanity checks
% FxH = GetVar(statefile, extrafile, {'VISCx_Um', 'VISCy_Um', ['Dx(1)/',num2str(dy*dz),'+Dy(2)/',num2str(dx*dz)]}, slice);
% FyH = GetVar(statefile, extrafile, {'VISCx_Vm', 'VISCy_Vm', ['Dx(1)/',num2str(dy*dz),'+Dy(2)/',num2str(dx*dz)]}, slice);
FxH = GetVar(statefile, extrafile, {'VISCx_Um', 'VISCy_Um', ['Dx(1)/',num2str(dy),'+Dy(2)/',num2str(dx)]}, slice);
FyH = GetVar(statefile, extrafile, {'VISCx_Vm', 'VISCy_Vm', ['Dx(1)/',num2str(dy),'+Dy(2)/',num2str(dx)]}, slice);
FxH = FxH./metric; %Variable z grid...
FyH = FyH./metric;

JFzH = bx.*FyH - by.*FxH;

%Numeric Friction Terms
FxN =  GetVar(statefile, diagfile, { 'Um_Advec', '(1)'}, slice)-UADVTERM;
FyN =  GetVar(statefile, diagfile, { 'Vm_Advec', '(1)'}, slice)-VADVTERM;

JFzN = bx.*FyN - by.*FxN;

% DIABATIC J Vectors
D = LHSb;

JBx = -OMEGAX.*D; JBy = -OMEGAY.*D; JBz = -OMEGAZ.*D;

% Horizontal Diabatic Terms
% DH = TtoB.*GetVar(statefile, extrafile, {'DFxE_TH', 'DFyE_TH', ['Dx(1)/',num2str(dy*dz),'+Dy(2)/',num2str(dx*dz)]}, slice);
DH = TtoB.*GetVar(statefile, extrafile, {'DFxE_TH', 'DFyE_TH', ['Dx(1)/',num2str(dy),'+Dy(2)/',num2str(dx)]}, slice);
DH = DH./metric;

JBzH = -OMEGAZ.*DH;

%Numeric Diabatic Terms
DN = -TtoB.*GetVar(statefile, diagfile, {'UDIAG1', '(1)'}, slice)-BADV;

JBzN = -OMEGAZ.*DN;


end