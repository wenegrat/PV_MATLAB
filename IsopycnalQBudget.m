%% 
% GENERATE ISOPYCNAL BUDGET
  clc;
 statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
% Parameters            
dx = 1000; dy = dx; dz = 5;
nx = 48; ny=nx; nz=100;
ts = 14400;
TtoB = 9.81.*2e-4;
tslice = [1 119];
slice={0, 0, 0, tslice};%100 120
sliceEta={0,0,[1 1],tslice};%251 271

isoT = [16.5 16.75];

%%
% CALCULATE TERMS
Q = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);
JBx = Q;
JBy = Q;
JBz = Q;
JFx = Q;
JFz = Q;
JFy = Q;
JAx = Q;
JAy = Q;
JAz = Q;
mask = Q;
for i=1:(tslice(end)-tslice(1)+1)
    disp(num2str(i))
   slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};
   [Q(:,:,:,i), JAx(:,:,:,i), JAy(:,:,:,i), JAz(:,:,:,i), ...
       JFx(:,:,:,i), JFy(:,:,:,i), JFz(:,:,:,i), JBx(:,:,:,i), JBy(:,:,:,i), JBz(:,:,:,i)] ...
       = calcQBudget(diagfile, statefile, etanfile, [nx, ny], slicetemp);
   THETA = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);
   mask(:,:,:,i) = (THETA>isoT(1)) & (THETA<isoT(2));
end

%%

