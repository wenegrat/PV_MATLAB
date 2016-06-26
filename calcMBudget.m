function [TEND, ADV, PRESS, DIFF, RES] = calcMBudget(statefile, diagfile, etanfile, divstr,dz, slice)
sliceEta = {slice{1}, slice{2}, [1 1], slice{4}};
TEND = GetVar(statefile, diagfile, {'TOTUTEND', '(1)/86400'}, slice);
ADV = GetVar(statefile, diagfile, {'Um_Advec', '(1)'}, slice);%Includes Coriolis
PRESS = GetVar(statefile, diagfile, {'Um_dPHdx', '(1)'}, slice);
dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.8*Dx(1)'},sliceEta));
[nx, ny, zL] = size(TEND);
PRESS = PRESS+permute(repmat((dEtadx), [1, 1,  zL]), [1 2 3]);
%Probably missing explicit flux terms here:
%Plus not sure about signs...
DIFF =  GetVar(statefile,diagfile, {'Um_Diss', 'VISrI_Um',['(1)-Dz(2)/',divstr]},slice);


RES = TEND-ADV-DIFF-PRESS;
end