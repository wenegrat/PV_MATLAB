function [out, dOutdt] =areaIntegrateJVecs(Jin, mask, dxdy, ts, vol)

dOutdt = squeeze(nansum(nansum(Jin(:,:,:).*mask(:,:,:)))).*dxdy;
out = cumtrapz(dOutdt).*ts./vol;

end