%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Theoretical Scalings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask = ones(nx, ny, nz, tslice(end)-tslice(1)+1);
vol = squeeze(sum(sum(sum(mask.*gridvol(:,:,1:nz,:)))));
gradb = TtoB.*(DPeriodic(THETA, dx, 'x').^2 + DPeriodic(THETA, dy, 'y').^2).^(1/2);

Hbl = GetVar(statefile, etanfile, {'KPPhbl', '(1)'}, sliceEta);

[nx, ny, nz, nt] = size(gradb);
Zfull = permute(repmat(Z, [1,nx, ny, nt]), [2 3 1 4]);

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


gradba = sqrt(squeeze(nanmean(nanmean(nanmean((gradb.^2).*masknan)))));

B0 = 9.81*2e-4*squeeze(Q0(1,1,1,1))./(1035*3994); %Only valid for constant Q0 in time.


%% FRICTIONAL THEORY SCALINGS

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
% DOMAIN AVERAGED FRICTION TERMS
ce = 0.08;
Jfeddy = +squeeze(ce*nanmean(nanmean(gradb(:,:,2,:).^2)).*nanmean(nanmean(Hbl)));

% [~, Jfdavg] = areaIntegrateJVecs(squeeze(Jfeddy), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);
% dOutdt = squeeze(nansum(nansum(Jin(:,:,:).*mask(:,:,:)))).*dxdy;
size(Jfeddy)
Jfdavg = squeeze(Jfeddy).*nx.*ny.*dx.*dy;

Jbdavg = 2*Jfdavg;

%% BUOYANCY THEORY SCALINGS

ce = 0.08;
JBsurf = -f0.*B0./Hbl;
JBeddy = 2*ce*gradb(:,:,2,:).^2.*Hbl;
JBtotp = JBsurf + JBeddy;


[~, Jbtotpa] = areaIntegrateJVecs(squeeze(JBtotp), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);
[~, Jbsurfa] = areaIntegrateJVecs(squeeze(JBsurf), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);
[~, Jbeddya] = areaIntegrateJVecs(squeeze(JBeddy), squeeze(masknan(:,:,2,:)), dx*dy, ts, vol);