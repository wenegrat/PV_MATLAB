function output = generateJVectorsRHS(pathpv,sliceT, OMEGAZ,OMEGAX, OMEGAY, Bx, By, alpha, beta, g, dx, dy,dz, mask, Q, Bz)

% U RHS
URHS = GetVarROMS(pathpv, 0, {'u_rhs', '(1)'}, sliceT);
% V RHS
VRHS = GetVarROMS(pathpv, 0, {'v_rhs', '(1)'}, sliceT);

% T RHS
TRHS = GetVarROMS(pathpv, 0, {'temp_rhs', '(1)'}, sliceT);

% S RHS
SRHS = GetVarROMS(pathpv, 0, {'salt_rhs', '(1)'}, sliceT);


JFZ = Bx.*VRHS - By.*URHS;
JFX = -Bz.*VRHS;
JFY = Bz.*URHS;

D   = g*alpha.*(TRHS)   + g*beta.*(SRHS);

% magb = sqrt(Bx.^2 + By.^2 + Bz.^2);
% FACT = (Q.*Bz./(magb.^2)); % See HM90 Eq. 4.1 and 4.2
JDZ = -(OMEGAZ).*D;
JDX = -OMEGAX.*D;
JDY = -OMEGAY.*D;

% XPROJ = Bx./magb; 
% YPROJ = By./magb;
% ZPROJ = Bz./magb;

output.JFZ = squeeze(JFZ(:,:,end-1));
output.JDZ = squeeze(JDZ(:,:,end-1));

output.Jfa = squeeze(nansum(nansum(squeeze(JFZ(:,:,end-1).*mask(:,:,end-1).*dx(:,:,end-1).*dy(:,:,end-1)))));
output.Jda = squeeze(nansum(nansum(squeeze(JDZ(:,:,end-1).*mask(:,:,end-1).*dx(:,:,end-1).*dy(:,:,end-1)))));

[nx ny nz] = size(D);
% X-Fluxes
xl = 2; xr = nx-1;
jfxl = squeeze(nansum(nansum(squeeze(JFX(xl,:,:).*mask(xl,:,:).*dy(xl,:,:).*dz(xl,:,:)))));
jfxr = squeeze(nansum(nansum(squeeze(JFX(xr,:,:).*mask(xr,:,:).*dy(xr,:,:).*dz(xr,:,:)))));
jdxl = squeeze(nansum(nansum(squeeze(JDX(xl,:,:).*mask(xl,:,:).*dy(xl,:,:).*dz(xl,:,:)))));
jdxr = squeeze(nansum(nansum(squeeze(JDX(xr,:,:).*mask(xr,:,:).*dy(xr,:,:).*dz(xr,:,:)))));

% Y-Fluxes
yf = 2; yb = ny-1;
jfyf = squeeze(nansum(nansum(squeeze(JFY(:,yf,:).*mask(:,yf,:).*dx(:,yf,:).*dz(:,yf,:)))));
jfyb = squeeze(nansum(nansum(squeeze(JFY(:,yb,:).*mask(:,yb,:).*dx(:,yb,:).*dz(:,yb,:)))));
jdyf = squeeze(nansum(nansum(squeeze(JDY(:,yf,:).*mask(:,yf,:).*dx(:,yf,:).*dz(:,yf,:)))));
jdyb = squeeze(nansum(nansum(squeeze(JDY(:,yb,:).*mask(:,yb,:).*dx(:,yb,:).*dz(:,yb,:)))));

 output.Jfa = output.Jfa + (jfxr - jfxl)+(jfyb-jfyf);
 output.Jda = output.Jda + (jdxr - jdxl)+(jdyb-jdyf);
end



