function z = compZ(path, type, zeta, theta_s, theta_b, hc, h)

% XXX - Guessing at these.
% Vtransform = 2; 
% Vstretching = 4; 



N = 50;

if ~type
    ts = 'r';
else
    ts = 'w';
end
% kgrid = 0; 
% column = 1;

% [nx ny] = size(h);
% z = NaN(nx, ny, N+kgrid);
% 
% for i=1:nx
% [z(i,:,:),s,C]=scoord(h, x, y, Vtransform, Vstretching, ...
%                         theta_s, theta_b, hc,             ...
%                         N, kgrid, column, i,0,0);
% end

% zeta = ncread(path, 'zeta');

z = zlevs(h, zeta, theta_s, theta_b, hc, N, ts, 2);
z = permute(z, [2 3 1]);
end