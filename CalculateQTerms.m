%% CALCULATE Q TERMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREALLOCATE FOR SPEED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);
% Qdir = Q;
JBx = Q;
JBy = Q;
JBz = Q;
JFx = Q;
JFz = Q;
JFy = Q;
JAx = Q;
JAy = Q;
JAz = Q;
JFzN = Q;
JBzN = Q;
JFzH = Q;
JBzH = Q;
THETA  = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);
OMEGAZs = NaN(nx, ny, tslice(end)-tslice(1) + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MULTIPLE PROCESSOR (speed increase factor of 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
tic
metric = permute(repmat(Z, [1, nx, ny, 1]), [2 3 1 4]);
for i=1:1:(tslice(end)-tslice(1)+1)
    disp(num2str(i));

   slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};
   [Q(:,:,:,i), JAx(:,:,:,i), JAy(:,:,:,i), JAz(:,:,:,i), ...
       JFx(:,:,:,i), JFy(:,:,:,i), JFz(:,:,:,i), JBx(:,:,:,i), JBy(:,:,:,i), JBz(:,:,:,i), ...
       JFzN(:,:,:,i), JBzN(:,:,:,i), JFzH(:,:,:,i), JBzH(:,:,:,i), OMEGAZs(:,:,i)]...
       = calcQBudgetD(diagfile, statefile, etanfile,extrafile, [nx, ny,nz, 1], slicetemp, dx, dy,dz);

   THETA(:,:,:,i) = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);
   
end
toc
