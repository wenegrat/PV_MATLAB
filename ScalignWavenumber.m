varw = squeeze(Bx(3:end-2,3:end-2,end-2,:));
varw = squeeze(JFT(3:end-2, 3:end-2,:));
[nx ny nt] = size(varw);
dx = 1500;
[psi, lambda] = sleptap(nx);  
ps = NaN(nt, ny, floor(nx/2)+1);
for j=1:nt
    disp(num2str(j));
for i=1:ny
    [k, ps(j, i,:)] = mspec(dx, squeeze(varw(:,i,j)), psi);
end
end

%%
psa = squeeze(nanmean(nanmean(ps)));
kn = k./(2*pi)*1000; %in 1/km
subplot(2,1,1)
semilogx(kn, kn.*psa, 'LineWidth', 2)
hold on;

% loglog(kn, (kn).^(-2)/1e12, '--');

hold off
xlabel('k (1/km)');
ylabel('PSD - J_{FT}');
title('Variance Preserving Spectra');
grid on

subplot(2,1,2)
semilogx(kn, cumtrapz(kn, psa)./trapz(kn,psa), 'LineWidth', 2);
grid on
xlabel('k (1/km)');
ylabel('Cumulative distribution');
hold on
yti = get(gca, 'YTick');
plot(1/3.*ones(size(yti)), yti);
hold off
set(gca, 'ylim', [0 1]);

set(gcf, 'Color' ,'w');