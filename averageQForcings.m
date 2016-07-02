function out =averageQForcings(input, mask, gridvol, ts)
[nx, ny, nz, nt] = size(input);
t = (1:nt).*ts;
vol = gridvol.*squeeze(sum(sum(sum(mask))));
% inputi = input.*mask;
% inputi = squeeze(nansum(nansum(nansum(inputi)))).*gridvol;
% inputi = cumtrapz(t,inputi);

inputt = cumtrapz(t, input, 4);
inputa = squeeze(nansum(nansum(nansum(inputt.*mask)))).*gridvol;
inputi = inputa;
out = inputi./vol;