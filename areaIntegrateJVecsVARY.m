function [out, dOutdt] =areaIntegrateJVecsVARY(Jin, mask, dxdy, ts, vol)
[nx, ny, nt] = size(Jin);
temp = NaN(nx, ny, nt);
for i=1:nx
    for j=1:ny
        temp(i,j,:) = cumtrapz(squeeze(Jin(i,j,:))).*ts;
    end
end
out = squeeze(nansum(nansum(temp.*mask))).*dxdy;
dOutdt = gradient(out, ts);
out = out-out(1);
out = out./vol;

    
    
end