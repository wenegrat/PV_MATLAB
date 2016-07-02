function [TEND, ADV, DIFF, RES] = calcBBudget(statefile, diagfile, etanfile, divstrh, divstrz,dz, slice, kppflag)
TtoB = 9.81.*2e-4;
TEND = GetVar(statefile, diagfile, {'TOTTTEND', '(1)/86400'}, slice);
ADV = GetVar(statefile, diagfile, {'UDIAG1', '(1)'}, slice);

%Using the built in flux diagnostics
% ADVx_TH = GetVar(statefile, diagfile, {'ADVx_TH', ['Dx(1)/',divstrz]}, slice);
% ADVy_TH = GetVar(statefile, diagfile, {'ADVy_TH', ['Dy(1)/',divstrz]}, slice);
% ADVr_TH = GetVar(statefile, diagfile, {'ADVr_TH', ['Dz(1)/',divstrh]}, slice);
% ADV= -ADVx_TH-ADVy_TH - ADVr_TH;

if kppflag
    DIFF = GetVar(statefile, diagfile, {'KPPg_TH','DFrI_TH', ['-Dz(1)/',divstr, '-Dz(2)/', divstr]}, slice);
else
    DIFF = GetVar(statefile, diagfile, {'DFxE_TH', 'DFyE_TH','DFrI_TH', ['-Dx(1)/',divstrz,'-Dy(2)/',divstrz,'-Dz(3)/', divstrh]}, slice);
end
TFLUX = GetVar(statefile, etanfile, {'TFLUX', '(1)'}, {0, 0, [1 1], slice{4}});
[nx ny nd] = size(DIFF);
TFLUXF = zeros(nx, ny, nd);
TFLUXF(:,:,1) = TFLUX(:,:,:);
Cw = 3994;		  
H = dz;
TFLUXF = TFLUXF./(1035*Cw*H);
DIFF = DIFF+TFLUXF;

RES = TEND-ADV-DIFF;
end