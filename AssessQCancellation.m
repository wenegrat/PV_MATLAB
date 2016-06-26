%%
% Calculate RHS terms in vorticity and buoyancy equatoin and check
% cancellation

% SEE VALLIS 4.58

function [residual, uterm, bterm] = AssessQCancellation(statefile, diagfile, slice)
bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);

%% 
% Calculate Absolute vorticity terms
% disp('Calculate Abs Vorticity Terms');
% disp('Vorticity Components');

OMEGAX = GetVar(statefile, diagfile, {'VVEL', ' - Dz(1)'}, slice);

OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);

OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL', '(1) + Dx(2) - Dy(3)'}, slice);

%%
%Velocity gradients
% disp('Velocity Terms')
ux = GetVar(statefile, diagfile, {'UVEL', 'Dx(1)'}, slice);
uy = GetVar(statefile, diagfile, {'UVEL', 'Dy(1)'}, slice);
uz = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);

vx = GetVar(statefile, diagfile, {'VVEL', 'Dx(1)'}, slice);
vy = GetVar(statefile, diagfile, {'VVEL', 'Dy(1)'}, slice);
vz = GetVar(statefile, diagfile, {'VVEL', 'Dz(1)'}, slice);

wx = GetVar(statefile, diagfile, {'WVEL', 'Dx(1)'}, slice);
wy = GetVar(statefile, diagfile, {'WVEL', 'Dy(1)'}, slice);
wz = GetVar(statefile, diagfile, {'WVEL', 'Dz(1)'}, slice);


xtermvel = bx.*(OMEGAX.*ux + OMEGAY.*uy + OMEGAZ.*uz);
ytermvel = by.*(OMEGAX.*vx + OMEGAY.*vy + OMEGAZ.*vz);
ztermvel = bz.*(OMEGAX.*wx + OMEGAY.*wy + OMEGAZ.*wz);

uterm = xtermvel + ytermvel + ztermvel;

%%
% disp('Buoyancy Terms')
u = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
v = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
w = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

[nx, ny, nz, nt] = size(u);
ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, nx, ny,nt ]), [2 3 1 4]);
% xgradbadv = Drv(dx, u.*bx, 'x');
% ygradbadv = Drv(dy, v.*by, 'y');
% zgradbadv = Drv(metric,w.*bz, 'z');

% xgradbadv = Drv(dx, u.*bx+v.*by+w.*bz, 'x');
% ygradbadv = Drv(dy,u.*bx+v.*by+w.*bz, 'y');
% zgradbadv = Drv(metric,u.*bx+v.*by+w.*bz, 'z');

xgradbadv = bx.*ux + by.*uy + bz.*uz;
ygradbadv = bx.*vx + by.*vy + bz.*vz;
zgradbadv = bx.*wx + by.*wy + bz.*wz;

bterm = OMEGAX.*xgradbadv + OMEGAY.*ygradbadv + OMEGAZ.*zgradbadv;
bterm = -bterm;
%%
% Found these to be O(e-29), therefore insignificant...

% Incompressible term
%  advdiv = ux +vy + wz;
%  icomp = bx.*OMEGAX.*advdiv + by.*OMEGAY.*advdiv + bz.*OMEGAZ.*advdiv;
%  % grad of curl term
%  omdiv = Drv(dx, OMEGAX, 'x') + Drv(dy, OMEGAY, 'y') + Drv(metric, OMEGAZ, 'z');
%  gcurl = bx.*u.*omdiv + by.*v.*omdiv + bz.*w.*omdiv;

%%
residual =  uterm+bterm;
% residual = residual.*mask;
end