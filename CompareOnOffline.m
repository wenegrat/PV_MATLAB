%%
% TEST ADV TERMS

%UVEL
%       DIRECT = GetVar(statefile, diagfile, { 'Um_Advec', '(1)'}, slice);
        U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
        Udir = ncread(diagfile, 'UVEL');
        Udir = squeeze(Udir(:,:,:,slice{4}(1):slice{4}(2)));
        %%
        plot(squeeze(U(indx, indy, indz, :)));
        hold on
        plot(squeeze(Udir(indx, indy, indz,:)));
        hold off
        %%
        V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
        W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);
        Ux = GetVar(statefile, diagfile, {'UVEL', 'Dx(1)'}, slice);
        Uy = GetVar(statefile, diagfile, {'UVEL', 'Dy(1)'}, slice);
        Uz = GetVar(statefile, diagfile, {'UVEL', 'Dz(1)'}, slice);
        OFFLINE = U.*Ux + V.*Uy + W.*Uz;
        
        ONLINE = GetVar(statefile, diagfile, {'UADVTERM', '(1)'}, slice);
        %%
        indx = 20; indy = 20; indz = 2;
        plot(squeeze(OFFLINE(indx,indy, indz,:)));
        hold on
                plot(squeeze(ONLINE(indx,indy, indz,:)));
        hold off
        %%
%VVEL

        Vx = GetVar(statefile, diagfile, {'VVEL', 'Dx(1)'}, slice);
        Vy = GetVar(statefile, diagfile, {'VVEL', 'Dy(1)'}, slice);
        Vz = GetVar(statefile, diagfile, {'VVEL', 'Dz(1)'}, slice);
        OFFLINE = U.*Vx + V.*Vy + W.*Vz;
        
        ONLINE = GetVar(statefile, diagfile, {'VADVTERM', '(1)'}, slice);
%         UDIAG5 = ncread(diagfile, 'UDIAG5');
%         UDIAG5 = squeeze(UDIAG5(:,:,:,(slice{4}(1):slice{4}(2))));
        
        %%
        indx = 20; indy = 20; indz = 2;
        plot(squeeze(OFFLINE(indx,indy, indz,:)), 'LineWidth', 2);
        hold on
                plot(squeeze(ONLINE(indx,indy, indz,:)));
%                 plot(squeeze(UDIAG5(indx, indy, indz,:)));
        hold off     
        
%%
%b
      
        bx = GetVar(statefile, diagfile, {'b', 'Dx(1)'}, slice);
        by = GetVar(statefile, diagfile, {'b', 'Dy(1)'}, slice);
        bz = GetVar(statefile, diagfile, {'b', 'Dz(1)'}, slice);
        OFFLINE = U.*bx + V.*by + W.*bz;
%         bdiag = GetVar(statefile, diagfile, {'UDIAG1', '(1)'}, slice);
        ONLINE = GetVar(statefile, diagfile, {'BADVTERM', '(1)'}, slice);
        %%
        indx = 10; indy = 2; indz = 200;
        figure
        plot(squeeze(OFFLINE(indx,indy, indz,:)));
        hold on
                plot(squeeze(ONLINE(indx,indy, indz,:)),'--');
                                plot(squeeze(bdiag(indx,indy, indz,:)),'--');

        hold off