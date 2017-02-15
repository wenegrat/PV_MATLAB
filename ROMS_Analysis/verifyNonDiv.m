mask = ones(size(Q));

xleft = xla(1) ; xright = xla(end) ;
dyl = repmat(squeeze(1./pn(xleft,:).'), [1 nz ]);
dyr = repmat(squeeze(1./pn(xright,:).'), [1 nz]);
dzl = squeeze(dz(xleft,:,:,end));
dzr = squeeze(dz(xright,:,:,end));
UL = sum(sum(squeeze(U(xleft,:,:)).*dyl.*dzl));
UR = sum(sum(squeeze(U(xright,:,:)).*dyr.*dzr));

deltaU = UR - UL;
% 
yfront = yla(1); yback = yla(end);
dxf = repmat(squeeze(1./pm(:,yfront)), [1 nz ]);
dxb = repmat(squeeze(1./pm(:,yback)), [1 nz ]);
dzf = squeeze(dz(:,yfront,:,end));
dzb = squeeze(dz(:,yback,:,end));
VF = sum(sum(squeeze(V(:,yfront,:)).*dxf.*dzf));
VB = sum(sum(squeeze(V(:,yback, :)).*dxb.*dzb));

deltaV = VB-VF;

WB = nansum(nansum(squeeze(W(:,:,2)).*1./pm.*1./pn));
deltaW = 0 - WB;

rem = deltaU + deltaV + deltaW

disp(['Percentage divergence: ', num2str(rem/vol(end))])