%% CALCULATE Q TERMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREALLOCATE FOR SPEED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);
Qdir = Q;
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
mask = Q;
FricDiv = Q;
AdvDiv = Q;
DiaDiv = Q;
THETA  = NaN(nx, ny, nz, tslice(end)-tslice(1)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINGLE PROCESSOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% incc =1;
% inc = incc-1;
% ztmp = ncread(statefile, 'Z');
% metric = permute(repmat(ztmp, [1, nx, ny, incc]), [2 3 1 4]);
% for i=1:incc:(tslice(end)-tslice(1)+1-inc)
%     disp(num2str(i));
% %     remaining = @(x)    remaining -1
% %     disp(['Remaining: ', num2str(remaining), '    (Processing: ', num2str(i),')'])
%    slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1+inc]};
%    [Q(:,:,:,i:i+inc), Qdir(:,:,:,i:i+inc), JAx(:,:,:,i:i+inc), JAy(:,:,:,i:i+inc), JAz(:,:,:,i:i+inc), ...
%        JFx(:,:,:,i:i+inc), JFy(:,:,:,i:i+inc), JFz(:,:,:,i:i+inc), JBx(:,:,:,i:i+inc), JBy(:,:,:,i:i+inc), JBz(:,:,:,i:i+inc), JFzN(:,:,:,i:i+inc), JBzN(:,:,:,i:i+inc)] ...
%        = calcQBudgetD(diagfile, statefile, etanfile, [nx, ny, incc], slicetemp, dx, dy);
% 
%    FricDiv(:,:,:,i:i+inc) = Drv(dx, JFx(:,:,:,i:i+inc), 'x') + Drv(dy, JFy(:,:,:,i:i+inc), 'y') + Drv(metric, JFz(:,:,:,i:i+inc), 'z');
% %    AdvDiv(:,:,:,i:i+inc) = Drv(dx, JAx(:,:,:,i:i+inc), 'x') + Drv(dy, JAy(:,:,:,i:i+inc), 'y') + Drv(metric, JAz(:,:,:,i:i+inc), 'z');
%    AdvDiv(:,:,:,i:i+inc) = DPeriodic(JAx(:,:,:,i:i+inc), dx, 'x') + DPeriodic(JAy(:,:,:,i:i+inc),dy, 'y') + Drv(metric, JAz(:,:,:,i:i+inc), 'z');
% 
%    DiaDiv(:,:,:,i:i+inc) = Drv(dx, JBx(:,:,:,i:i+inc), 'x') + Drv(dy, JBy(:,:,:,i:i+inc), 'y') + Drv(metric, JBz(:,:,:,i:i+inc), 'z');
%    THETA(:,:,:,i:i+inc) = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);
% %    [resid(:,:,:,i:i+inc), uterm(:,:,:,i:i+inc), bterm(:,:,:,i:i+inc)] = AssessQCancellation(statefile, diagfile, slicetemp);
% end

%Do end bit
% i = i +incc;
% slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(2)]};
% metric = permute(repmat(ztmp, [1, nx, ny, slicetemp{4}(2)-slicetemp{4}(1)+1]), [2 3 1 4]);
% 
%    [Q(:,:,:,i:end), Qdir(:,:,:,i:end), JAx(:,:,:,i:end), JAy(:,:,:,i:end), JAz(:,:,:,i:end), ...
%        JFx(:,:,:,i:end), JFy(:,:,:,i:end), JFz(:,:,:,i:end), JBx(:,:,:,i:end), JBy(:,:,:,i:end), JBz(:,:,:,i:end), JFzN(:,:,:,i:end), JBzN(:,:,:,i:end)] ...
%        = calcQBudgetD(diagfile, statefile, etanfile, [nx, ny, slicetemp{4}(2)-slicetemp{4}(1)+1], slicetemp, dx, dy);
% 
%    FricDiv(:,:,:,i:end) = Drv(dx, JFx(:,:,:,i:end), 'x') + Drv(dy, JFy(:,:,:,i:end), 'y') + Drv(metric, JFz(:,:,:,i:end), 'z');
% %    AdvDiv(:,:,:,i:i+inc) = Drv(dx, JAx(:,:,:,i:i+inc), 'x') + Drv(dy, JAy(:,:,:,i:i+inc), 'y') + Drv(metric, JAz(:,:,:,i:i+inc), 'z');
%    AdvDiv(:,:,:,i:end) = DPeriodic(JAx(:,:,:,i:end), dx, 'x') + DPeriodic(JAy(:,:,:,i:end),dy, 'y') + Drv(metric, JAz(:,:,:,i:end), 'z');
% 
%    DiaDiv(:,:,:,i:end) = Drv(dx, JBx(:,:,:,i:end), 'x') + Drv(dy, JBy(:,:,:,i:end), 'y') + Drv(metric, JBz(:,:,:,i:end), 'z');
%    THETA(:,:,:,i:end) = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MULTIPLE PROCESSOR (speed increase factor of 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
metric = permute(repmat(Z, [1, nx, ny, 1]), [2 3 1 4]);
parfor i=1:1:(tslice(end)-tslice(1)+1)
    disp(num2str(i));

   slicetemp = {slice{1}, slice{2}, slice{3}, [tslice(1)+i-1 tslice(1)+i-1]};
   [Q(:,:,:,i), Qdir(:,:,:,i), JAx(:,:,:,i), JAy(:,:,:,i), JAz(:,:,:,i), ...
       JFx(:,:,:,i), JFy(:,:,:,i), JFz(:,:,:,i), JBx(:,:,:,i), JBy(:,:,:,i), JBz(:,:,:,i), JFzN(:,:,:,i), JBzN(:,:,:,i)] ...
       = calcQBudgetD(diagfile, statefile, etanfile, [nx, ny, 1], slicetemp, dx, dy);

   FricDiv(:,:,:,i) = DPeriodic(JFx(:,:,:,i), dx, 'x') + DPeriodic(JFy(:,:,:,i),dy, 'y') + Drv(metric, JFz(:,:,:,i), 'z');
   AdvDiv(:,:,:,i) = DPeriodic(JAx(:,:,:,i), dx, 'x') + DPeriodic(JAy(:,:,:,i),dy, 'y') + Drv(metric, JAz(:,:,:,i), 'z');
   DiaDiv(:,:,:,i) = DPeriodic(JBx(:,:,:,i), dx, 'x') + DPeriodic(JBy(:,:,:,i),dy, 'y') + Drv(metric, JBz(:,:,:,i), 'z');

   THETA(:,:,:,i) = GetVar(statefile, diagfile, {'THETA', '(1)'}, slicetemp);
end