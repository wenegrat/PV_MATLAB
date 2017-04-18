function [DPDY, VADV, VCORI, VTEND, IFRICV, TY, AX, AY, AZ, TYY]= returnMomTerms(diagfile, statefile, etanfile,extrafile, sizes, slice, dx, dy,dz )
nx = sizes(1); ny = sizes(2); nz=sizes(3); nt=sizes(4);
sliceEta = {slice{1}, slice{2}, [1 1], slice{4}};

% Load Velocities
U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

ztmp = ncread(statefile, 'Z');
ztmp = ztmp(1:nz);
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);
zL = length(ztmp);


dpdy = GetVar(statefile,diagfile, {'Um_dPHdx','(1)'},slice);
dEtady = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.81*(1)'},sliceEta));  
dEtady = DPeriodic(dEtady, dx, 'x');
DPDY = permute(repmat((dEtady), [1, 1, 1, zL]), [1 2 4 3]) + dpdy;

Vx = DPeriodic(U, dx, 'x');
Vy = DPeriodic(U, dy, 'y');
Vz = Drv(metric, U, 'z');
VADV= -U.*Vx - V.*Vy - W.*Vz ;

VCORI =  GetVar(statefile, diagfile, {'Um_Cori', '(1)'}, slice);

VTEND = GetVar(statefile, diagfile, {'TOTUTEND', '(1)/86400'}, slice);

IFRICV = VTEND - VCORI - DPDY - VADV;

% VADV = -W.*Vz;% XX-TESTING THIS

T = GetVar(statefile, diagfile, {'THETA', '(1)'}, slice);
TY = DPeriodic(T, dy, 'y');
TYY = DPeriodic(TY, dy, 'y');
AX = -U.*Vx; AY = -V.*Vy; AZ = -W.*Vz;

end


        
        
