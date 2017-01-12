function [out, dOutdt] =areaIntegrateJVecsROMS(Jin, mask, dx, dy, ts, vol)

dOutdt = squeeze(nansum(nansum(Jin(:,:,:).*mask(:,:,:).*dx.*dy)));
out = cumtrapz(dOutdt).*ts./vol;

end