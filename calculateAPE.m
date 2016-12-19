function [APE PE SC] = calculateAPE(THETA,X, Y, Z,T)
[nx, ny, nz, nt] = size(THETA);

x = linspace(X(1), X(end), 100);
y = linspace(Y(1), Y(end), 100);
z = linspace(Z(1), Z(end), 100);

THETAred = interpn(X, Y, Z, T, THETA, x, y, z, T);
[nx, ny, nz, nt] = size(THETAred);

% First sort temp profiles
TS = NaN(size(THETAred));
for i=1:nt;
%     thets = interp3(X, Y, Z, squeeze(THETA(:,:,:,i)), x, y, z);
    tvec = reshape(THETAred(:,:,:,i), nx*ny*nz,1);
    tvecs = sort(tvec, 'descend');
    TS(:,:,:,i) = reshape(tvecs, nx, ny, nz);
end   
% [X1 Y1 Z1 T1] = ndgrid(X, Y, Z, T);
% TSi = interpn(x, y, z, T, TS, X1, Y1, Z1, T1);

Tref = 16;
rho = 1035.*(1 - 2e-4.*(THETAred-Tref));
% zfull = repmat(z, [100 1 100 nt]);
% PE = trapz(z, 9.81.*rho.*zfull, 3);
rhom = squeeze(nansum(nansum(rho)));
PE = trapz(z, 9.81.*repmat(z.', [1 nt]).*rhom);

rhos = 1035.*(1-2e-4.*(TS - Tref));
% APE = PE - trapz(z, 9.81.*rhos.*zfull, 3);
rhosm = squeeze(nansum(nansum(rhos)));
APE = PE - trapz(z, 9.81.*repmat(z.', [1 nt]).*rhosm);

b0 = 9.81*2e-4*squeeze(100)./(1035*3994);
for j=1:nt
zstarind = find(rhom(1,j) >= rhosm(:,j), 1, 'last');
zstar = z(zstarind+1);
SC(j) = 1035.*b0.*zstar;
end

% for i =1
% SurfChange = rho.*100.*
% APE = s`queeze(nansum(nansum((APE))));
% PE = squeeze(nansum(nansum((PE))));
end