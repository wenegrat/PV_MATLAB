%%
dx = 1000; dy = 1000; dz = 2.5;

% Following JAS article
bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);
b = GetVar(statefile, diagfile, {'b', '(1)'},slice);
%%
% Equation 1b
% Is my vorticity incompressible?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OMEGAX = GetVar(statefile, diagfile, {'WVEL', 'VVEL', 'Dy(1) - Dz(2)'}, slice);
OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'WVEL', 'Dz(1) - Dx(2)'}, slice);
% OMEGAX = GetVar(statefile, diagfile, {'VVEL', ' - Dz(1)'}, slice);
% OMEGAY = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);
OMEGAZ = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL', '(1) + Dx(2) - Dy(3)'}, slice);
ztmp = ncread(statefile, 'Z');
[nx, ny, nz, nt] = size(OMEGAX);
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);
div = Drv(dx, b.*OMEGAX, 'x') + Drv(dy, b.*OMEGAY, 'y') + Drv(metric, b.*OMEGAZ, 'z');
Q = bx.*OMEGAX + by.*OMEGAY + bz.*OMEGAZ;

res1b = (Q - div)./Q;
% Find residual is not always small locally, but seems small in a volume
% integrated budget. Normalized mean is O(-3e-3). Although some instances
% of larger values...
figure('Name', 'Eq. 1b residual')
plot(squeeze(nanmean(nanmean(nanmean(res1b)))))

%%
% Equation 6a
% Can I balance my budget like this?
% Which Q is balanced here (direct or normal)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
