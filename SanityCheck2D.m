statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc'; extrafile = 'extra.nc';
kppfile = 'kppdiags.nc';

TtoB = 9.81.*2e-4;
f0 = 1e-4;

b = TtoB.*outputFull.T;
bx = DPeriodic(b, 500,'x');
by = DPeriodic(b, 500, 'y');

[nx ny nz nt] = size(b);

ztmp = ncread(statefile, 'Z');
ztmp = ztmp(1:nz);
metric = permute(repmat(ztmp, [1, nx, ny, nt]), [2 3 1 4]);
bz = Drv(metric, b, 'z');
Z = ncread(statefile, 'Z');
Zl = ncread(statefile, 'Zl');
magb = abs(bx+1i.*by);

Hkpp = GetVar(statefile, etanfile, {'KPPhbl', '(1)'}, {0, 0, [1 1], 0});
W = GetVar(statefile, diagfile, {'WVEL' ,'(1)'}, {0, 0, 0, 0});
V = GetVar(statefile, diagfile, {'VVEL' ,'(1)'}, {0, 0, 0, 0});
Vz = Drv(metric, V(:,:,:,1:end-1), 'z');
U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, {0, 0, 0, 0});
Uz = Drv(metric, U(:,:,:,1:end-1), 'z');
Uza = Uz - (-by./f0);
Bo = 9.81*2e-4*squeeze(25)./(1035*3994);
difft = GetVar(statefile, kppfile, {'KPPdiffT', '(1)'}, {0, 0, 0, 0});
diffm = GetVar(statefile, kppfile, {'KPPviscA', '(1)'}, {0, 0, 0, 0});
%%
Rib = bz.*f0.^2./(by.^2);
fzone = 1:200;
zl= 1:50;
% vprime = squeeze(V(2,fzone, zl, :)) - squeeze(repmat((nanmean(nanmean(V(2,fzone, zl, :)))), [1 length(fzone) length(zl) 1]));
% wprime = squeeze(W(2,fzone, zl, :)) - squeeze(repmat((nanmean(nanmean(W(2,fzone, zl, :)))), [1 length(fzone) length(zl) 1]));
% bprime = squeeze(b(2,fzone, zl, :)) - squeeze(repmat((nanmean(nanmean(b(2,fzone, zl, :)))), [1 length(fzone) length(zl) 1]));

%Average in y only
vprime = squeeze(V(2,fzone, zl, :)) - squeeze(repmat((nanmean(V(2,fzone, zl, :))), [1 length(fzone) 1 1]));
wprime = squeeze(W(2,fzone, zl, :)) - squeeze(repmat((nanmean(W(2,fzone, zl, :))), [1 length(fzone) 1 1]));
bprime = squeeze(b(2,fzone, zl, :)) - squeeze(repmat((nanmean(b(2,fzone, zl, :))), [1 length(fzone) 1 1]));

vpwp = squeeze(nanmean((vprime.*wprime)));
bpwp = squeeze(nanmean((bprime.*wprime(:,:,1:end-1))));
vpbp = squeeze(nanmean(bprime.*vprime(:,:,1:end-1)));
GSP = -vpwp(:,1:end-1).*squeeze(nanmean((by(2,fzone,zl,:)./f0)));

% bpwpt = squeeze
plot(nanmean(GSP))
hold on
plot(nanmean(bpwp));
hold off
%%
vm = max(diffm, [], 3);
% vm = squeeze(nanmean(vm(2,fzone,1,:)));
EBFg = 1./f0.^2.*squeeze(nanmean(nanmean((diffm(2,fzone,1:16,1:end-1).*by(2,fzone,1:16,:).^2))));
tl = 300:400;
plot(nanmean(GSP(:,tl), 2), Zl)
hold on
plot(nanmean(bpwp(:,tl), 2), Zl)
plot(nanmean(bpwp(:,tl), 2)+nanmean(GSP(:,tl), 2), Zl, '--')
xt = get(gca, 'XTick');
plot(xt, squeeze(nanmean(nanmean(Hkpp(2,fzone,1,tl)))).*-1*ones(size(xt)))
yt = get(gca, 'Ytick');
plot(nanmean(EBFg(tl)).*ones(size(yt)), yt);
hold off
%%

    
ms = double(permute(repmat(outputFull.Z, [1 nx ny nt]), [2 3 1 4]) > -repmat(Hkpp(:,:,:,2:end), [1 1 nz 1]));

Bzh = bz.*ms;
bzh = NaN(nx,ny, 1, nt);
diffh = bzh;
byh = bzh;
B = findfirst(~ms, 3); %% XX-Should check that this is working correctly.
B = B-2;
B(B>nz) = nz;
        for xi=1:nx;
            for yi =1:ny;
                for tt = 1:nt
                if B(xi, yi, nt)<nz
                bzh(xi,yi,1,tt) = bz(xi,yi, B(xi,yi,tt), tt);
                diffh(xi,yi,1,tt) = difft(xi,yi, B(xi,yi,tt), tt);
                byh(xi,yi,1,tt) = by(xi, yi, B(xi,yi,tt), tt);
                end
                end
            end
        end



dbzmld = diffh.*bzh;
%%
fzone = 1:120;
% fzone = 40:60;
fzone = 1:200;
fzone = 75:120;
jfm = squeeze(nanmean(outputFull.JFz(2,fzone,2,:)));
jft = squeeze(nanmean(Hkpp(2,fzone,1,2:end).*magb(2,fzone,2,:).^2));
jbm = squeeze(nanmean(outputFull.JBz(2,fzone,2,:)));
jbt = squeeze(nanmean(Hkpp(2,fzone,1,2:end).*magb(2,fzone,2,:).^2));
jbs = squeeze(nanmean(f0.*Bo./Hkpp(2,fzone, 1,2:end)));

effb = -squeeze(-0.03.*jft.*squeeze(nanmean(Hkpp(2,fzone,1,2:end)))./f0);
jbs = squeeze(nanmean(f0.*(Bo./Hkpp(2,fzone, 1, 2:end) + 1.*dbzmld(2,fzone,1,:)./Hkpp(2,fzone, 1, 2:end))));
jbs = 1.2.*squeeze(nanmean(f0.*(Bo./Hkpp(2,fzone,1,2:end))));
plot(-jfm);
hold on
plot(jft);
hold off

tl = 1:length(jfm);
cf = regress(jfm(tl), jft(tl))
% regress(output.dJf, output.dJfea)

% cds = regress(jbm(tl), jbs(tl))
cd =regress(jbm(tl)-1.*jbs(tl), jbt(tl))
% regress(output.dJb, output.dJbea)

cf + cd
%%
plot(jbm)
hold on
plot(jbs)
hold off

%%
corr(jbm(tl), jbs+cd.*jbt)
bd = regress(jbm(tl), jbs+cd.*jbt)
b = regress(jbm(tl), [jbs(tl) jbt(tl)])

corr(jbm(tl), b(1).*jbs +  b(2).*jbt)
%%
subplot(2,1,1)
plot(squeeze(nanmean(outputFull.JFz(2,:,2,:),4)))
subplot(2,1,2)
plot(squeeze(nanmean(dbzmld(2,:,1,:), 4))./Bo)
% plot(squeeze(nanmean(bzh(2,:,1,:), 4))./(8.*f0).^2)
% plot(squeeze(nanmean(bzh(2,:,1,:), 4)).*f0.^2./squeeze(nanmean(byh(2,:,1,:),4)).^2);
% plot(squeeze(nanmean(bzh(2,:,1,:), 4)).^2.*f0.^2./squeeze(nanmean(byh(2,:,1,:),4)).^2);

% plot(squeeze(nanmean(wb(2,:,1,:), 4))./Bo)

% plot(squeeze(nanmean(diffh(2,:,1,:), 4)))

%%
Havg = squeeze(nanmean(Hkpp(2,fzone,1,ti)));

fzone = 100:120;
wb =V(:,:,:,1:end-1).*b;%(b-repmat(nanmean(b(:,fzone,:,:),2), [1 ny 1 1]));
% wb = W(:,:,:,1:end-1).*bz(:,:,:,1:end);
ti = 100;
plot(squeeze(nanmean(wb(2,fzone,:,ti)))./(Bo), outputFull.Z)
hold on;
plot(squeeze(nanmean(difft(2,fzone,:,ti).*bz(2,fzone,:,ti)))./Bo, outputFull.Z);
xt = get(gca, 'XTick');
plot(xt, -ones(size(xt)).*Havg);
plot(1, 0, 'x')
hold off
%%
Havg = squeeze(nanmean(Hkpp(2,fzone,1,ti)));

fzone =75:125;
vb =V(:,:,:,1:end-1).*by;%(b-repmat(nanmean(b(:,fzone,:,:),2), [1 ny 1 1]));
wb = W(:,:,:,1:end-1).*bz(:,:,:,1:end);
ebf = .06.*Havg.*(Bo*Havg).^(1/3)./f0.*squeeze(nanmean(nanmean(by(2,:,:,ti).*ms(2,:,:,ti))))./Havg;
ti = 100;
plot(squeeze(nanmean(vb(2,fzone,:,ti))), outputFull.Z)
hold on;
plot(squeeze(nanmean(wb(2,fzone,:,ti))), outputFull.Z)

% plot(squeeze(nanmean(difft(2,fzone,:,ti).*bz(2,fzone,:,ti)))./Bo, outputFull.Z);
xt = get(gca, 'XTick');
plot(xt, -ones(size(xt)).*Havg);
yt = get(gca, 'YTick');
plot(ebf.*ones(size(yt)), yt);
% plot(1, 0, 'x')
hold off
grid on
%%
plot(squeeze(difft(2,100,:,100)), outputFull.Z);
hold on
plot(squeeze(diffm(2,100,:,100)), outputFull.Z);
hold off