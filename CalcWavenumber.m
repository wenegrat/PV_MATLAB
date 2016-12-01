%%
% Calculate KE
slice={0, 0, [1 5], tslice};
% KE = GetVar(statefile, diagfile, {'UVEL', 'VVEL', '(1).^2 + (2).^2'}, slice);
U = GetVar(statefile, diagfile, {'UVEL' '(1).^2'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
[nx, ny, nz, nt] = size(U);
% Up = U - repmat(nanmean(U,4), [1 1 1 nt]);
% Vp = V - repmat(nanmean(V,4), [1 1 1 nt]);
Up = U - repmat(nanmean(U,1), [nx 1 1 1]);
Vp = V - repmat(nanmean(V,1), [nx 1 1 1]);
KE = Up.^2 + Vp.^2;
%%
% Calculate wavenumber spectra


[psi, lambda] = sleptap(nx);    
ps = NaN(nt, ny, floor(nx/2)+1);
for j=1:nt
for i=1:ny
    [k, ps(j, i,:)] = mspec(dx, squeeze(KE(:,i,1,j)), psi);
end
end

%%
figure
%Plots
kn = k./(2*pi)*1000;
% loglog(k, ps);
loglog(kn, (kn).^(-2)/1e1, '--');

hold on
for i=(120:60:480)
loglog(kn, squeeze(nanmean(ps(i,:,:),2)), 'LineWidth', 2);
% loglog(kn, squeeze(ps(i,end,:)), 'LineWidth', 2);

end
hold off
hold on;

loglog(kn, (kn).^(-3)/1e2, '--');
yti = get(gca, 'YTick');
loglog(1/1.5.*ones(size(yti)), yti);
hold off
grid on
legend('Model', 'k^{-2}', 'k^{-3}');    
xlabel('km^{-1}');
ylabel('KE PSD');
set(gcf, 'Color', 'w');
%%
% pcolor(kn, 1:nt, log10(squeeze(nanmean(ps,2))));
% shading interp