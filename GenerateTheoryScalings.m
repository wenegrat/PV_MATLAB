%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Theoretical Scalings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tslice = [1 1339];
% slice={0, 0, [1 5], tslice};
% sliceEta={0,0,[1 1],tslice};
mask = ones(nx, ny, nz, tslice(end)-tslice(1)+1);
vol = squeeze(sum(sum(sum(mask.*gridvol(:,:,1:nz,:)))));

% Q0 = ncread(etanfile, 'TFLUX', slice);
% Q0 = GetVar(statefile, etanfile, {'TFLUX', '(1)'}, sliceEta);
% THETA = GetVar(statefile, diagfile, {'THETA', '(1)'}, slice);

gradb = TtoB.*(DPeriodic(THETA, dx, 'x').^2 + DPeriodic(THETA, dy, 'y').^2).^(1/2);

% ZETA = GetVar(statefile, extrafile, {'momVort3', '(1)'}, sliceEta);
% ZETA = GetVar(statefile, diagfile, {'VVEL', 'UVEL', 'Dx(1)-Dy(2)'}, slice);
%%
Hbl = GetVar(statefile, etanfile, {'KPPhbl', '(1)'}, sliceEta);
% ghat = GetVar(statefile, kppfile, {'KPPg_TH',['(1)/',num2str(dx*dy)]}, slice);
% bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
% by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
% bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);
%%
[nx, ny, nz, nt] = size(gradb);
% magbgrad = bx.^2 + by.^2;
% K = GetVar(statefile, kppfile, {'KPPdiffT', '(1)'}, slice);
Zfull = permute(repmat(Z, [1,nx, ny, nt]), [2 3 1 4]);

%%

masknan = Zfull(:,:,1:nz,:)>-repmat(Hbl, [1, 1, nz, 1]);

masknan = double(masknan);
masknan(masknan<1) = NaN;

%Try adding a mask to the frontal zone only.
% [nx ny nz nt] = size(by);
% bym = nanmedian(reshape(abs(by(:,:,1,:)), nx*ny, nt));
% bymm = repmat(bym, [nx 1 ny nz]);
% bymm = permute(bymm, [1 3 4 2]);
% fzonemask = abs(by) > 1.1.*bymm;
% masknan = masknan.*fzonemask;

% gradba = squeeze(nanmean(nanmean(nanmean(gradb.*masknan))));
% gradban = gradba./gradba(1); % ONLY APPROXIMATE
gradba = sqrt(squeeze(nanmean(nanmean(nanmean((gradb.^2).*masknan)))));
% gradbvar = (squeeze(nanvar(nanvar(nanvar((gradb.*masknan))))));
% gradbmed = nanmedian(reshape(gradb.*masknan, [nx*ny*nz, nt]));

% Hbla = squeeze(nanmean(nanmean(Hbl)));
B0 = 9.81*2e-4*squeeze(Q0(1,1,1,1))./(1035*3994);

% rho = 1035*(1-2e-4.*(THETA(:,:,1,:)-16.5));
% B0 = 9.81*2e-4*squeeze(Q0(1,1,1,1))./(rho*3994);

% JfT = Hbla.^(4/3).*abs(B0).^(1/3).*gradba.*2;
% JfT = Hbla.^(1/3).*abs(B0).^(1/3).*gradbvar;

%% FRICTIONAL THEORY SCALINGS
% Jftot = -squeeze(gradb(:,:,1,:)).*squeeze(gradb(:,:,1,:)).*squeeze(Hbl).^(4/3).*abs(squeeze(B0)).^(1/3)/2;
% Jftot = -squeeze(gradb(:,:,1,:)).*squeeze(gradb(:,:,1,:)).*sqrt(squeeze(Hbl).^(4/3).*abs(squeeze(B0)).^(1/3)/(f0*2));

%Spatially varying
Vg = Hbl.*gradb(:,:,1,:)./f0;
wstar = (abs(B0)*Hbl).^(1/3);
% nukpp = .1*wstar.*Hbl;

%Note that max(G) ~ 4/27 at sigma = 1/3
% so max(G)*0.4 ~ 0.06
Vttw = .061*wstar.*Vg./(f0.*Hbl);
% Vttw = .1*wstar.*Vg./(f0*Hbl);
Jfgeo =-squeeze(f0*Vttw.*gradb(:,:,2,:));

ce = 0.08;
Jfeddy = +squeeze(ce*gradb(:,:,2,:).^2.*Hbl);

% INCLUDE EDDY TERMS?
Jftot = Jfgeo + Jfeddy;

% Jftot = -sqrt(0.03).*sqrt(Hbl.*wstar./f0).*gradb(:,:,2,:).^2;

[~, Jftota] = areaIntegrateJVecs(squeeze(Jftot), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);
[~, Jfeddya] = areaIntegrateJVecs(squeeze(Jfeddy), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);
[~, Jfgeoa] = areaIntegrateJVecs(squeeze(Jfgeo), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);


%%
% Spatially Averaged
% Vg = squeeze(nanmean(nanmean(Hbl.*masknan(:,:,2,:)))).*squeeze(nanmean(nanmean(gradb(:,:,2,:).*masknan(:,:,2,:))))./f0;
% wstar = squeeze(nanmean(nanmean((abs(B0)*Hbl.*masknan(:,:,2,:)).^(1/3))));
% % nukpp = .1*wstar.*Hbl;
% 
% %Note that max(G) ~ 4/27 at sigma = 1/3
% % so max(G)*0.4 ~ 0.06
% Vttw = .061*wstar.*Vg./(f0.*Hbla);
% Jftota =- squeeze(f0*Vttw.*squeeze(nanmean(nanmean(gradb(:,:,2,:)))));
% 
% [~, Jftota] = areaIntegrateJVecs(-squeeze(f0*permute(repmat(Vttw, [1, nx, ny, 1]), [2 3 4 1]).*gradb(:,:,2,:)), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);
% Jftot = nukpp.*gradb(:,:,2,:).^2./(f0.*Hbl);
% Hbl0 = permute(repmat(squeeze(nanmean(nanmean(Hbl))), [1 nx ny 1]), [2 3 4 1]);
% wstar = (abs(B0)*Hbl0).^(1/3);
% nukpp = abs(.1.*wstar.*Hbl0);
% de = sqrt(2.*nukpp./f0);
% Jftot = -de.*gradb(:,:,2,:).^2./2;

% Jftota = squeeze(nansum(nansum(Jftot.*squeeze(masknan(:,:,2,:)*dx*dy))));
%% BUOYANCY THEORY SCALINGS
% JBtot = B0.*(f0+ZETA)./Hbl; % Based on bulk arguments over boundary layer
% JBtotp = B0.*f0./(Hbl); %Corr = 0.73 % Bulk arguments but not including zeta
% JBtotp = (f0+ZETA).*(B0-K(:,:,2,:).*bz(:,:,2,:)+TtoB.*ghat(:,:,2,:)); % All 3 terms corr = 0.95
% JBtotp = (f0+ZETA).*(B0-K(:,:,2,:).*bz(:,:,2,:)); %Only Surf and Local corr = 0.89
ce = 0.08;
% ce = 0.06;
% mfk = sqrt(by.^2 - bx.^2) ;
% mfk = nanmean(gradb(:,:,1:4,:),3);

% Ek = .1*wstar./(f0*Hbl.^2);
JBsurf = -f0.*B0./Hbl;
JBeddy = 2*ce*gradb(:,:,2,:).^2.*Hbl;

% JBtotp = -f0.*(B0 - 2*ce.*mfk(:,:,:,:).^2.*Hbl.^2./f0)./Hbl; %Surf and <w'b'> Fox-Kemper corr = 0.92.
% JBtotp = -f0.*(B0*0 - 2*ce.*mfk(:,:,:,:).^2.*Hbl.^2./f0)./Hbl; %Surf and <w'b'> Fox-Kemper corr = 0.92.
JBtotp = JBsurf + JBeddy;

% kh = squeeze(nanmean(nanmean(ce.*bz(:,:,2,:).*Hbl.^2./f0)));
% JBtotp = f0.*(ce.*mfk(:,:,1,:).*bz(:,:,2,:).*Hbl.^2./f0)./dy; %Surf and <w'b'> Fox-Kemper corr = 0.92.
% JBtotp = -f0.*(mfk.*f
% JBtotp = -f0.*(B0 - 0.5*mfk(:,:,:,:).*bz(:,:,2,:).*Hbl.^2./f0)./Hbl;
%B0~m^2/s^3,   
% disp('here')
%  JBtotp = ce.*f0.*(magbgrad(:,:,2,:).*Hbl.^2./f0)./Hbl; % Only <w'b'> corr = 0.95
% JBtotp = -(f0).*B0./Hbl;
% JBtotp = 0.06.*JBtotp; % FK;
% Ri = f0.^2.*nanmean(bz(:,:,:,:).*masknan, 3)./magbgrad(:,:,2,:);
% mu = (Hbl.*abs(B0)).^(1/3)./(f0*Hbl);
% JBtotp = f0.*(magbgrad(:,:,2,:).*Hbl.^2./f0).*(1-mu.^2).^(-2)./Hbl ;

% Jbtotpa = squeeze(nanmean(nanmean(JBtotp.*masknan(:,:,2,:)))).*vol;
[~, Jbtotpa] = areaIntegrateJVecs(squeeze(JBtotp), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);
[~, Jbsurfa] = areaIntegrateJVecs(squeeze(JBsurf), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);
[~, Jbeddya] = areaIntegrateJVecs(squeeze(JBeddy), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);