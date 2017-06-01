function output= generateJVectors(pathts, pathuv, sliceT,Bx, By,Bz, rho0, g, OMEGAX, OMEGAY, OMEGAZ, mask, dxf, dyf,dz, alpha, beta)
dx = squeeze(dxf(:,:,1)); dy = squeeze(dyf(:,:,1));

thmix = GetVarROMS(pathts, 0, {'temp_hmix', '(1)'}, sliceT);
tvmix = GetVarROMS(pathts, 0, {'temp_vmix', '(1)'}, sliceT);
shmix = GetVarROMS(pathts, 0, {'salt_hmix', '(1)'}, sliceT);
svmix = GetVarROMS(pathts, 0, {'salt_vmix', '(1)'}, sliceT);
sforc = GetVarROMS(pathts, 0, {'salt_forc', '(1)'}, sliceT);
tforc = GetVarROMS(pathts, 0, {'temp_forc', '(1)'}, sliceT);

[nx, ny, nz] = size(thmix);
alpha = repmat(alpha, [1 1 nz]); % NOTE THIS IS NOT EXACT, but much faster...
beta = repmat(beta, [1 1 nz]);
D     = g*alpha.*(thmix + tvmix+tforc)   + g*beta.*(shmix + svmix+sforc);
output.Jd = -OMEGAZ.*D;
JBx = -OMEGAX.*D; JBy = -OMEGAY.*D;

uhmix = GetVarROMS(pathuv, 0, {'u_hmix', '(1)'}, sliceT);
uhdif = GetVarROMS(pathuv, 0, {'u_hdiff', '(1)'}, sliceT);
uvmix = GetVarROMS(pathuv, 0, {'u_vmix', '(1)'}, sliceT);
vhmix = GetVarROMS(pathuv, 0, {'v_hmix', '(1)'}, sliceT);
vvmix = GetVarROMS(pathuv, 0, {'v_vmix', '(1)'}, sliceT);
vhdif = GetVarROMS(pathuv, 0, {'v_hdiff', '(1)'}, sliceT);

Fx = uhmix + uvmix+uhdif;
Fy = vhmix + vvmix +vhdif;
JFx = -Bz.*Fy; JFy = Bz.*Fx;

output.Jf = Bx.*Fy - By.*Fx;

% Vertical Terms
output.Jfa = squeeze(nansum(nansum(output.Jf(:,:,end-1).*squeeze(mask(:,:,end-1)).*dx.*dy)));
output.Jda = squeeze(nansum(nansum(output.Jd(:,:,end-1).*squeeze(mask(:,:,end-1)).*dx.*dy)));

% Zonal Terms
dJFx_l = squeeze(nansum(nansum(squeeze(JFx(2,:,:).*mask(2,:,:)).*squeeze(dyf(2,:,:).*dz(2,:,:)))));
dJFx_r = squeeze(nansum(nansum(squeeze(JFx(end-1,:,:).*mask(end-1,:,:)).*squeeze(dyf(end-1,:,:).*dz(end-1,:,:)))));
output.dJFxA = dJFx_r - dJFx_l;
dJBx_l = squeeze(nansum(nansum(squeeze(JBx(2,:,:).*mask(2,:,:)).*squeeze(dyf(2,:,:).*dz(2,:,:)))));
dJBx_r = squeeze(nansum(nansum(squeeze(JBx(end-1,:,:).*mask(end-1,:,:)).*squeeze(dyf(end-1,:,:).*dz(end-1,:,:)))));
output.dJBxA = dJBx_r - dJBx_l;
% Meridional Terms
dJFy_f = squeeze(nansum(nansum(squeeze(JFy(:,2,:).*mask(:,2,:)).*squeeze(dxf(:,2,:).*dz(:,2,:)))));
dJFy_b = squeeze(nansum(nansum(squeeze(JFy(:,end-1,:).*mask(:,end-1,:)).*squeeze(dxf(:,end-1,:).*dz(:,end-1,:)))));
output.dJFyA = dJFy_b - dJFy_f;
dJBy_f = squeeze(nansum(nansum(squeeze(JBy(:,2,:).*mask(:,2,:)).*squeeze(dxf(:,2,:).*dz(:,2,:)))));
dJBy_b = squeeze(nansum(nansum(squeeze(JBy(:,end-1,:).*mask(:,end-1,:)).*squeeze(dxf(:,end-1,:).*dz(:,end-1,:)))));
output.dJByA = dJBy_b - dJBy_f;
end