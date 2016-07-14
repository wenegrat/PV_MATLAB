% nx = 96; ny = 192;
% dz = 3;

nx = 80; ny=100;
dz = 3;
nz = 100;
 tslice = [1 1300];
slice={0, 0, 0, tslice};%100 120
divstrh = '122500';
divstrz = '1050';
kpp = true;
TEND = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);
ADV = TEND;
DIFF = TEND;
RES = TEND;
TENDZ = TEND;
DIFFZ = TEND;
ADVZ = TEND;
RESZ = TEND;
ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, nx, ny, 1]), [2 3 1 4]);
bz = TEND;
for i = 1:(tslice(end)-tslice(1)+1);
    disp(num2str(i));
    
   slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};

    [TEND(:,:,:,i), ADV(:,:,:,i), DIFF(:,:,:,i), RES(:,:,:,i)] = calcBBudget(statefile, diagfile, etanfile,kppfile,extrafile, divstrh, divstrz,dz, slicetemp, kpp);
%     TENDZ(:,:,:,i) = Drv(metric, squeeze(TEND(:,:,:,i)), 'z');
%     ADVZ(:,:,:,i) = Drv(metric, squeeze(ADV(:,:,:,i)), 'z');
%     DIFFZ(:,:,:,i) = Drv(metric, squeeze(DIFF(:,:,:,i)), 'z');
%     RESZ(:,:,:,i) = Drv(metric, squeeze(RES(:,:,:,i)), 'z');
%     bz(:,:,:,i) = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slicetemp);
end
RHS = ADV+DIFF;

%%
ind = 198;
tmean = nanmean(nanmean(TEND(:,:,ind)));

figure
pcolor(squeeze(RES(:,:,ind)./tmean));
colorbar

%%
indx = 1; indy = 45; indz = 1;

figure
plot(squeeze(TEND(indx, indy, indz,:)));
hold on
plot(squeeze(ADV(indx, indy, indz,:)));
plot(squeeze(DIFF(indx, indy, indz,:)));
plot(squeeze(RES(indx, indy, indz,:)), 'LineWidth', 2);
plot(squeeze(RHS(indx, indy, indz,:)), 'LineWidth', 2, 'LineStyle', '--');
% plot(squeeze(abgt(indx, indy, indz,:)));
hold off
legend('TEND', 'ADV', 'DIFF', 'RES', 'RHS')
%%
% plot(squeeze(TENDZ(indx, indy, indz,:)));
% hold on
% plot(gradient(squeeze(bz(indx, indy, indz, :))./(9.81.*2e-4), 1800));
% hold off
% 
% %%
% plot(cumtrapz(squeeze(TENDZ(indx, indy, indz,:))).*1800);