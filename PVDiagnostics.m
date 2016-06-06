%%

pv = GetVar('~/MITgcm/projects/PVINJECT3/state_0000000000.nc', 0, {'PV', '(1)'},{0, 0, 0, 0});
pvv = GetVar('~/MITgcm/projects/PVINJECT3/state_0000000000.nc', 0, {'PVv', '(1)'},{0, 0, 0, 0});

jbz = XXX;
jbf = XXX;

%%
time = ncread('~/MITgcm/projects/PVINJECT3/state_0000000000.nc', 'T');
%%
PV = squeeze(nansum(nansum(nansum(pv))));
PVV= squeeze(nansum(nansum(nansum(pvv))));

plot(time./(86400), PV);
hold on
plot(time./86400, PVV, '--');
hold off