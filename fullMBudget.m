nx = 48; ny = 96; nz = 200;
dz = 2.5;
tslice = [40 110];
slice={0, 0, 0, tslice};%100 120
divstr = '1e6';
TEND = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);
ADV = TEND;
DIFF = TEND;
RES = TEND;
PRESS = TEND;
ztmp = ncread(statefile, 'Z');
metric = permute(repmat(ztmp, [1, nx, ny, 1]), [2 3 1 4]);
bz = TEND;
for i = 1:(tslice(end)-tslice(1)+1);
    disp(num2str(i));
    
   slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};

    [TEND(:,:,:,i), ADV(:,:,:,i), PRESS(:,:,:,i), DIFF(:,:,:,i), RES(:,:,:,i)] = calcMBudget(statefile, diagfile, etanfile, divstr,dz, slicetemp);
  
end
RHS = ADV+PRESS+DIFF;
%%
ind = 2;
tmean = nanmean(nanmean(TEND(:,:,ind)));

figure
pcolor(squeeze(RES(:,:,ind)));
colorbar

%%
indx = 20; indy = 48; indz = 1;

figure
plot(squeeze(TEND(indx, indy, indz,:)));
hold on
plot(squeeze(ADV(indx, indy, indz,:)));
plot(squeeze(PRESS(indx, indy, indz,:)));
plot(squeeze(DIFF(indx, indy, indz,:)));
plot(squeeze(RES(indx, indy, indz,:)), 'LineWidth', 2);
plot(squeeze(RHS(indx, indy, indz,:)), 'LineWidth', 2, 'LineStyle', '--');

hold off
legend('TEND', 'ADV', 'PRESS','DIFF', 'RES')
%%
plot(squeeze(TENDZ(indx, indy, indz,:)));
hold on
plot(gradient(squeeze(bz(indx, indy, indz, :))./(9.81.*2e-4), 1800));
hold off

%%
plot(cumtrapz(squeeze(TENDZ(indx, indy, indz,:))).*1800);

%%
scatter(squeeze(ADV(indx, indy, indz, :)), squeeze(PRESS(indx, indy, indz,:)))