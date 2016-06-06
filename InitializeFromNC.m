%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Generate Relevant Diagnostic Variables

%Following Thomas 2005
%Convert from b = -alpha*g*Temp = 9.81*2e-4*Temp
statefile = 'state.nc';
diagfile = 'diag.nc';
etanfile = 'etan.nc';
disp('Jbz')
slice = {0, 0, 0, 0};
D = 9.81.*2e-4.*GetVar(statefile, diagfile, {'TOTTTEND', 'UDIAG1', '(1)/86400 - (2)'}, slice);

%OMEGA_A = fk + Grad x u = i*[ dyW - dzV] + j*[dzU-dxW] + k*[f + dxV -dyU]
%q = f*OMEGA_A dot Grad B = f* [ (dyW-dzV)dxB + (dzU-dxW)dyB + (f + dxV-dyU)*dzB

% qt = (f + dxV - dyU)dzB -dzV*dxB + dzU*dyB


Omegaz = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL', '(1)+Dx(2)-Dy(3)'},slice);
Jbz = -1e-4.*Omegaz.*D; %Vertical component of Buoyancy J Vector

Omegay = GetVar(statefile, diagfile, {'UVEL', 'WVEL', 'Dz(1)-Dx(2)'}, slice);
Jby = -1e-4.*Omegay.*D; %Meridional component of Buoyancy J Vector

% Jbz = GetVar(statefile, diagfile, {'f_1e-4', 'VVEL', 'UVEL','TOTTTEND','UDIAG1', ...
%              '-(1).*((1) + Dx(2)-Dy(3)).*9.81.*2e-4.*((4)/86400-(5))'},{0, 0, 0, 0});
[nx ny nz nt ] = size(Jbz);
%%
disp('bx')   
bx = GetVar(statefile, diagfile, {'b',...
      ' Dx(1)  '},slice);
disp('Fy')   
Fy = GetVar(statefile, diagfile, {'TOTVTEND','Vm_dPHdy','Vm_Advec',...
      ' (1)/86400 - (2) -(3) '},slice);
disp('etay')   
dEtady = squeeze(GetVar(statefile,etanfile, {'ETAN',...
      '-9.8*Dy(1)'},{0, 0,[1 1], slice{4}}));
Fy = Fy - permute(repmat(squeeze(dEtady),[1, 1, 1, nz]), [1 2 4 3]);
disp('by')
by = GetVar(statefile, diagfile, {'b',...
      ' Dy(1)  '},slice);
disp('Fx')
Fx = GetVar(statefile, diagfile, {'TOTUTEND','Um_dPHdx','Um_Advec',...
      ' (1)/86400 - (2) -(3) '},slice);
disp('etax')
dEtadx = squeeze(GetVar(statefile, etanfile, {'ETAN',...
      '-9.8*Dx(1)'},{0, 0,[1 1], slice{4}}));
Fx = Fx - permute(repmat(squeeze(dEtadx),[1, 1, 1, nz]), [1 2 4 3]);
Jfz = 1e-4*(bx.*Fy - by.*Fx); %Vertical component of friction J vector

bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);
%Assuming vertical component of F vector is zero from hydrostatic balance.
Jfy = 1e-4.*(bz.*Fx); % Meridional component of friction J vector
%%
% Calculate PV
q = GetVar(statefile, diagfile, {'PV',...
      '(1)'},slice);
qv = GetVar(statefile, diagfile,{'PVv',...
      '(1)'},slice);
  
qfull = GetVar(statefile, diagfile, {'f_1e-4', 'WVEL', 'VVEL', 'b', 'UVEL', '(1).*( (Dy(2)-Dx(3)).*Dx(4) + (Dz(5)-Dx(2)).*Dy(4) + ((1) + Dx(3)-Dy(5)).*Dz(4))'}, slice);
%%
UVEL = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
VVEL = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
WVEL = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);

% {'UVEL' 'VVEL' 'f_1e-4' 'b' ['((3)+Dx(2)-Dy(1)).*Dz(4)-Dz(2).*' ...
%                         'Dx(4)+Dz(1).*Dy(4)']};
%%
% Try pointwise budget
indj = 40; indk = 20;
qvar = qfull;1e-4.*qv;
% qt = gradient( squeeze(5*nansum(squeeze(qfull(indj, indk, :,:)))), ts); %single column, integrated w.r.t z
qt = squeeze(5*nansum(squeeze(qvar(indj, indk, :,:))))- squeeze(5*nansum(squeeze(qvar(indj, indk, :,1))));
jtemp =  squeeze(-Jfz(indj, indk,:,:) - Jbz(indj, indk,:,:));
jz = jtemp(1,:) - jtemp(end,:);
jz = cumtrapz(jz)*ts;
plot(qt);
hold on
plot(-jz);
hold off

%%
% Jbza = -cumtrapz(squeeze(nansum(nansum(nansum(Jbz))))).*3600;
% Jfza = -cumtrapz(squeeze(nansum(nansum(nansum(Jfz))))).*3600;
Lx=48e3;    
Ly = 48e3;
Lz = 500;
ts = 3600;
area = Lx*Ly;
domx = 1:nx; domy = 1:ny;

%Averaging Horizontally
Jbah = squeeze(nansum(nansum(Jbz(domx,domy,:,:))))*1000*1000./area;
Jfah = squeeze(nansum(nansum(Jfz(domx,domy,:,:))))*1000*1000./area;

Jbzt = Jbz;
Jbzt(~isfinite(Jbzt)) = 0;%not really correct...
Jbah = squeeze(trapz(1000*(1:nx), Jbzt)./Lx);
Jbah = squeeze(trapz(1000*(1:ny), Jbah)./Ly);

Jfzt = Jfz;
Jfzt(~isfinite(Jfzt)) = 0;
Jfah = squeeze(trapz(1000*(1:nx), Jfzt)./Lx);
Jfah = squeeze(trapz(1000*(1:ny), Jfah)./Ly);

Jfyah = 5.*1000/(Lx*Ly).*squeeze(nansum(nansum(Jfy(domx,:,:,:), 1), 3));
Jbyah = 5.*1000/(Lx*Ly).*squeeze(nansum(nansum(Jby(domx, :,:,:), 1), 3));
%JAyah = 5.*1000/(Lx*Ly).*squeeze(nansum(nansum(VVEL(domx, :,:,:).*qvar(domx,:,:,:), 1), 3));

%Integrating the vertical difference in time.
Jbza = cumtrapz(Jbah(1,:) - Jbah(end,:)).*ts;
Jfza = cumtrapz(Jfah(1,:) - Jfah(end,:)).*ts;

%For meridional component of frictional flux integrate merid diff in time.
Jfya = -cumtrapz(Jfyah(domy(end),:) - Jfyah(domy(1),:)).*ts;
Jbya = -cumtrapz(Jbyah(domy(end),:) - Jbyah(domy(1),:)).*ts;
%Jaya = -cumtrapz(JAyah(domy(end),:) - JAyah(domy(1),:)).*ts;
%%
% PV
%Pick a PV Variable
qvar = 1e-4*q;
%Averaging horizontally
qvar(~isfinite(qvar)) = 0;
qvarah = squeeze(trapz(1000*(1:nx), qvar))./Lx;
qvarah = squeeze(trapz(1000*(1:ny), qvarah))./Ly;

% qvarah = squeeze(nansum(nansum(qvar(domx,domy,:,:))))*1000*1000./area;

%Differencing in time (ie define a DeltaQ.
qvarahdelta = qvarah(:,:) - repmat(qvarah(:,1), 1, length(qvarah));
%Integrate vertically from top to bottom (maybe needs to be negative?)
% qa = 1e-4.*nansum(qvarahdelta(2:end,:))*5;
qa = trapz(5*-1*(0:99), qvarahdelta);
%%

plot(qa); hold on; plot(-Jbza, '--k'); plot(-Jfza); plot(-(Jbza + Jfza), '--');plot(-Jfya-Jbya); hold off
legend('Q_v', 'JBZA', 'J_F', 'Sum', 'Jy_F', 'Jy_B');

%%
scatter(qa, 2*(Jbza+Jfza))
grid on