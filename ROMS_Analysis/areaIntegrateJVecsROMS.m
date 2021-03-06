function [out, dOutdt] =areaIntegrateJVecsROMS(Jin, mask, dx, dy, ts, vol)
% disp('j')
dOutdt = squeeze(nansum(nansum(Jin(:,:,:).*mask(:,:,:).*dx.*dy)));
% dOutdt = smooth(dOutdt, 5);
out = cumtrapz(dOutdt).*ts./vol;


% out = squeeze(nansum(nansum(cumsum(Jin(:,:,:),3).*ts.*mask(:,:,:).*dx.*dy)))./vol;
end