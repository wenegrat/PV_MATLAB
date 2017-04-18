

[psi, lambda] = sleptap(nt);    
variab = squeeze(output.Qt);
% variab = squeeze(nanmean(nanmean(Vz(2,fzone,1:20,:))));
variab = squeeze(nanmean(bpwp(1:16,1:end)).');
% variab = squeeze(nanmean(nanmean(Vz(2,100:105, 2:3,:).*by(2,100:105,2:3,:))));
[sig, ps] = mspec(2, variab, psi);
sign = sig./(2*pi);

feff = nanmean(sqrt(f0.*squeeze(nanmean(nanmean(nanmean((outputFull.Q(:,:,16:28,:))))))./squeeze(nanmean(nanmean(nanmean((bz(:,:,16:28,:))))))));
figure
semilogx(sign, sig.*ps)
yt = get(gca, 'YTick');
hold on; plot(f0./(2*pi).*3600.*ones(size(yt)), yt); hold off
hold on; plot(f0./(2*pi).*3600.*ones(size(yt))./2, yt); hold off
hold on; plot(f0./(2*pi).*3600.*ones(size(yt)).*2, yt); hold off

hold on; plot(abs(feff)./(2*pi).*3600.*ones(size(yt)), yt); hold off