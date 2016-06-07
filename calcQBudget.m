function [Q, JADVx, JADVy, JADVz, JFx, JFy, JFz, JBx, JBy, JBz] = calcQBudget(diagfile, statefile, etanfile, sizes, slice )
nx = sizes(1); ny = sizes(2);
sliceEta = {slice{1}, slice{2}, [1 1], slice{4}};
TtoB = 9.81.*2e-4;

% disp(slice{1}); disp(slice{2}); disp(slice{3}); disp(slice{4});
ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, nx, ny, 1]), [2 3 1 4]);
zL = length(ztmp);

dpdx = GetVar(statefile,diagfile, {'Um_dPHdx','(1)'},slice);
dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.8*Dx(1)'},sliceEta));
% if ndims(dEtadx)==2
%     disp('NDIMS == 2');
%     temp = dEtadx;
%     clear dEtadx;
%     dEtadx(1,:,:) = zeros(size(temp));
% end
dpdx = permute(repmat((dEtadx), [1, 1, 1, zL]), [1 2 4 3]) + dpdx;
dpdy = GetVar(statefile,diagfile, {'Vm_dPHdy','(1)'},slice);
dEtady = squeeze(GetVar(statefile, etanfile, {'ETAN','-9.8*Dy(1)'},sliceEta));  
% if ndims(dEtady)==2`
%     disp('NDIMS == 2');
%     temp = dEtady;
%     clear dEtady;
%     dEtady(1,:,:) = temp;
% end
dpdy = permute(repmat((dEtady), [1, 1, 1, zL]), [1 2 4 3]) + dpdy;

%%
% Calculate buoyancy gradients
% disp('Calculate B gradients')
bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
% if ndims(squeeze(bx))==3; bx = zeros(size(bx)); end
by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);

%% 
% Calculate Absolute vorticity terms
% disp('Calculate Abs Vorticity Terms');

%ignore w derivatives for consistency with hydrostatic approx.
OMEGAX = GetVar(statefile, diagfile, {'WVEL', 'VVEL', 'Dy(1) - Dz(2)'}, slice);
OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'WVEL', 'Dz(1) - Dx(2)'}, slice);

% OMEGAX = GetVar(statefile, diagfile, { 'VVEL', ' - Dz(1)'}, slice);
% OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'Dz(1) '}, slice);
OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL', '(1) + Dx(2) - Dy(3)'}, slice);
% if ndims(squeeze(OMEGAZ))==3; 
%     OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'UVEL', '(1)  - Dy(2)'}, slice);
%     OMEGAY = Drv(metric, GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice),'z');
% end
%%
Q = OMEGAX.*bx + OMEGAY.*by + OMEGAZ.*bz; %A more direct definition of Q.

% Calculate J Vectors
% disp('Calculate J Vectors');
U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);
% 
% %ADVECTIVE TERMS
JADVx = U.*Q; 
JADVz = W.*Q;
JADVy = V.*Q;
%FRICTION TERMS
Fx = GetVar(statefile, diagfile, {'TOTUTEND','Um_Advec', ' (1)/86400 - (2)'},slice);
Fx = Fx - dpdx;
Fy = GetVar(statefile, diagfile, {'TOTVTEND', 'Vm_Advec', '(1)/86400 - (2)'}, slice);
Fy = Fy - dpdy;

JFx = -bz.*Fy; JFy = bz.*Fx; JFz = bx.*Fy - by.*Fx;

%DIABATIC TERMS
D = TtoB*GetVar(statefile, diagfile, {'TOTTTEND', 'UDIAG1', '(1)/86400 - (2)'}, slice);
JBx = -OMEGAX.*D; JBy = -OMEGAY.*D; JBz = -OMEGAZ.*D;
%%

end