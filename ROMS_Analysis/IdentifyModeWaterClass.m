

%%
xpos = 10;
figure
for i=100:nt
    
    contour(1:ny, zm, squeeze(Rho(xpos, :,:,i)).'); shading interp
    colorbar
    set(gca, 'ylim', [-1500 0])
    pause()
    
end


%%
figure
xi = 10; yi = 20; ti=120;
yl = [-800 0];

subplot(1,3,1)
plot(squeeze(Rho(xi,yi,:,ti)), zm)
set(gca, 'ylim', yl);
grid on

subplot(1,3,2)
plot(squeeze(Tf(xi,yi,:,ti)), zm)
set(gca, 'ylim', yl);
grid on

subplot(1,3,3)
plot(squeeze(Bz(xi,yi,:,ti)), zm)
set(gca, 'ylim', yl);
grid on

%%
rhovec  = reshape(rho, nx*ny*nz*nt,1);
zvec = reshape(dz, nx*ny*nz*nt,1);
rhovec(rhovec<1000) = NaN;
zvec = zvec(isfinite(rhovec));
rhovec = rhovec(isfinite(rhovec));

edges=1020:.1:1040;
edges=linspace(min(rhovec), max(rhovec), 100);
[trash bin]=histc(rhovec,edges);
count=accumarray(bin,zvec(:));

subplot(1,2,1)
bar(edges,count)
%%
tvec  = reshape(Tf, nx*ny*nz*nt,1);

edges=linspace(min(tvec),max(tvec),100);
[trash bin]=histc(tvec,edges);
count=accumarray(bin,zvec(:));

subplot(1,2,2)
bar(edges,count)