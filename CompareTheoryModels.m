% foldernames(1,:) = 'GS_025_6F';
% foldernames(2,:) = 'GS_200_6F';
% foldernames(3,:) = 'GS_100_3F';
% foldernames(4,:) = 'GS_200_3F';
% foldernames(5,:) = 'GS_025_2F';%
% foldernames(6,:) = 'GS_200_2F';%
% foldernames(7,:) = 'GS_025_4F';%
% foldernames(8,:) = 'GS_200_4F';
cd /data/thomas/jacob13/PARAMSPACE/
%% FIRST LOAD ALL OUTPUT FLAT FILES
clear foldernames
foldernames(1,:) = 'GS_025_1F';
foldernames(2,:) = 'GS_025_2F';
foldernames(3,:) = 'GS_025_4F';
foldernames(4,:) = 'GS_025_6F';
foldernames(5,:) = 'GS_100_1F';%
foldernames(6,:) = 'GS_100_2F';%
foldernames(7,:) = 'GS_100_4F';%
foldernames(8,:) = 'GS_100_6F';%
foldernames(9,:) = 'GS_200_1F';%
foldernames(10,:) = 'GS_200_2F';%
foldernames(11,:) = 'GS_200_4F';
foldernames(12,:) = 'GS_200_6F';
% foldernames(10,:) = 'GS_025_1F';
% foldernames(11,:) = 'GS_100_1F';

[nr, ~] = size(foldernames);

jbsc = 1.2; %Empirical coefficient from 1D case.

% clear foldernames
% foldernames(1,:) = 'GS_100_4F_LEK';
% nr = 1;

Jfm = NaN(nr, 1000);
Jbm = Jfm; Jft = Jfm;  Jbt = Jfm;
Jea = Jfm; Jga=Jfm;
Jbs = Jfm;
Jbe = Jfm;
Jbe_zeta =Jfm;
Jbs_zeta = Jfm;
Jbm_zeta = Jfm;
for i=1:nr;
    filename = ['./',foldernames(i,:),'/', foldernames(i,:),'_OutputsFlat.mat'];
    load(filename);
    nt = length(output.dJf);
    Jfm(i,1:nt) = output.dJf;
    Jbm(i,1:nt) = output.dJb;
    Jbm_zeta(i,1:nt) = output.dJb_zeta;
%     Jbm(i,1:12) = NaN;
%     Jfm(i,1:12) = NaN;
    Jft(i,1:nt) = output.dJfea;
%     Jea(i,1:nt) = -output.dJfea;
%     Jga(i,1:nt) = output.dJfga; % This is the TTW scaling.

%     Jbt(i,1:nt) = output.dJbsa+output.dJbea;
    Jbs(i,1:nt) = output.dJbsa*jbsc;
    Jbe(i,1:nt) = output.dJbea;
    Jbs_zeta(i,1:nt) = output.dJbsa_zeta.*jbsc;
    Jbe_zeta(i,1:nt) = output.dJbea_zeta;
    
    cs = strsplit(foldernames(i,:), '_');
    legstring(i,:) = ['Q_o: -', cs{2}, ', M^2: (', cs{3},')^2'];
end
%% FIND BEST FIT COEFFICIENT
[nruns nt] = size(Jfm);
mvec = reshape(Jfm, nruns*nt,1);
tvec = reshape(abs(Jft), nruns*nt, 1);
% gvec = reshape(abs(Jga), nruns*nt, 1);
mask = isfinite(mvec);
Fco = abs(regress(mvec(mask), tvec(mask)));
disp(['Best fit Frictional Coefficient: ', num2str(Fco)]);
mvec = reshape(Jbm-Jbs, nruns*nt,1);
svec = reshape(abs(Jbe), nruns*nt, 1);
mask = isfinite(mvec);
Dco = regress(mvec(mask), svec(mask));
disp(['Best fit Diabatic Coefficient: ', num2str(Dco)]);

%% FIND BEST FIT COEFFICIENT ZETA
[nruns nt] = size(Jfm);

mvec = reshape(Jbm_zeta-Jbs, nruns*nt,1);
svec = reshape(abs(Jbe), nruns*nt, 1);
mask = isfinite(mvec);
DcoZ = regress(mvec(mask), svec(mask));
disp(['Best fit Diabatic Zeta Coefficient: ', num2str(DcoZ)]);
%%
% disp('-------------------- Bootstrap Version ----------------------')
% mvec = reshape(Jfm, nruns*nt,1);
% tvec = reshape(abs(Jft), nruns*nt, 1);
% bootdat = [mvec(mask),tvec(mask)];
% bootb = bootstrp(1000, @(x) regress(x(:,1),x(:,2:end)),bootdat);
% bootb = sort(bootb);
% Fco = abs(median(bootb));
% Ll = bootb(25); Ul = bootb(975);
% disp(['Best fit Frictional Coefficient: ', num2str(Fco)]);
% disp(['CI Limits Frictional Coefficient: ', num2str(Ll),' - ', num2str(Ul)]);
% mvec = reshape(Jbm-Jbs, nruns*nt,1);
% svec = reshape(abs(Jbe), nruns*nt, 1);
% mask = isfinite(mvec);
% % Dco = regress(mvec(mask), svec(mask));
% bootdat = [mvec(mask),svec(mask)];
% bootb = bootstrp(1000, @(x) regress(x(:,1),x(:,2:end)),bootdat);
% bootb = sort(bootb);
% Dco = median(bootb);
% Ll = bootb(25); Ul = bootb(975);
% disp(['Best fit Diabatic Coefficient: ', num2str(Dco)]);
% disp(['CI Limits Diabatic Coefficient: ', num2str(Ll),' - ', num2str(Ul)]);
%%
[nruns nt] = size(Jfm);
xl = [1e-6 1e-1];
% H = 150; f0 = 1e-4;Q = 100; gradb = (3*f0).^2;
% B = 9.81*2e-4*Q./(1035*3994);
% Vg = H.*gradb./f0;
% wstar = (abs(B)*H).^(1/3);
% Vttw = .061*wstar.*Vg./(f0.*H);
% Jfi = f0*Vttw*gradb*dx*dy*nx*ny;

meanfac = 1/12;
% meanfac =1;
figure
subplot(1,2,1)
hold on
map = colormap(parula(8));
for i=1:nruns;
     if mod(i,4)==1 %1F
       set(gca, 'ColorOrderIndex', 1);

     elseif mod(i,4) ==2 %2F
                 set(gca, 'ColorOrderIndex', 2);

%          mark = 'o'
     elseif mod(i,4) ==3 %4F
                 set(gca, 'ColorOrderIndex', 3);
     else    %6F
                 set(gca, 'ColorOrderIndex', 4);

%          mark = 's';
     end
     if i<5
         mark = 'd';
     elseif i>4 && i<9
          mark = 'o';
     else
         mark = 's';
     end
%     if i==4; set(gca, 'ColorOrderIndex', 5); end %Skip Purple

    % Note that I am getting rid of the ce factors here, and using the
    % regression coefficent found at end of this script. Better approach
    % would be to save flat files with no ce coefficient, and then find the
    % regression coefficient here before plotting.
    scatter(abs(imresize(Jfm(i,:), meanfac)),Fco.*abs(imresize(Jft(i,:), meanfac)),mark,  'filled', 'MarkerEdgeColor', 'k');
%     scatter(abs(Jfm(i,1:1/meanfac:end)),Fco.*abs(Jft(i,1:1/meanfac:end)),mark,  'filled', 'MarkerEdgeColor', 'k');

end
%     s =scatter(Jfi, Jfi, 'x', 'LineWidth', 3);
set(gca, 'xscale', 'log', 'yscale', 'log');
xt = get(gca, 'XTick');
plot(xt, xt, 'k');
plot(xt, 2*xt,'--k')
plot(xt, 0.5*xt, '--k')
set(gca, 'xlim', xl, 'ylim', xl);
hold off
mvec = reshape(Jfm, nruns*nt,1);
tvec = reshape(-Jft, nruns*nt, 1);
mask = isfinite(mvec+tvec);
cr = corr(mvec(mask), tvec(mask));
title(['Frictional PV Flux,     Corr: ', num2str(cr,2)]);
grid on
% xlabel('Model  $(m^3s^{-4})$'); ylabel('Scaling $(m^3s^{-4})$');
xlabel('$|J_F|$ $(m^3s^{-4})$'); ylabel('$c_{F} H |\nabla_H b|^2$ $(m^3s^{-4})$');

set(gca, 'FontSize', 16);
subplot(1,2,2)
hold on
for i=1:nruns;
      if mod(i,4)==1 %1F
       set(gca, 'ColorOrderIndex', 1);

     elseif mod(i,4) ==2 %2F
                 set(gca, 'ColorOrderIndex', 2);

%          mark = 'o'
     elseif mod(i,4) ==3 %4F
                 set(gca, 'ColorOrderIndex', 3);
     else    %6F
                 set(gca, 'ColorOrderIndex', 4);

%          mark = 's';
     end
     if i<5
         mark = 'd';
     elseif i>4 && i<9
          mark = 'o';
     else
         mark = 's';
     end
    scatter(abs(imresize(Jbm(i,:), meanfac)), abs(imresize(Jbs(i,:)+Dco.*Jbe(i,:), meanfac)), mark, 'filled', 'MarkerEdgeColor', 'k');
%         scatter(abs(Jbm(i,1:1/meanfac:end)),abs(Jbs(i,1:1/meanfac:end)+Dco.*Jbe(i,1:1/meanfac:end)),mark,  'filled', 'MarkerEdgeColor', 'k');

end
set(gca, 'xscale', 'log', 'yscale', 'log');
set(gca, 'xlim', xl, 'ylim', xl);

xt = linspace(xl(1),xl(end), 100);
plot(xt, xt, 'k');
plot(xt, 2*xt,'--k')
plot(xt, 0.5*xt, '--k')

hold off
mvec = reshape(Jbm, nruns*nt,1);
tvec = reshape(Jbs+Dco.*Jbe, nruns*nt, 1);
mask = isfinite(mvec+tvec);
cr = corr(mvec(mask), tvec(mask));
title(['Diabatic PV Flux,     Corr: ', num2str(cr,2)]);
grid on
% xlabel('Model $(m^3s^{-4})$'); ylabel('Scaling $(m^3s^{-4})$');
xlabel('$|J_D|$ $(m^3s^{-4})$'); ylabel('$c_sf B_o/H + c_{D} H |\nabla_H b|^2$ $(m^3s^{-4})$');

%legend(legstring(1,:), legstring(2,:),legstring(3,:),legstring(4,:), legstring(5,:), legstring(6,:),legstring(7,:),legstring(8,:), legstring(9,:), 'location', 'SouthEastOutside');
% l = legend('$\nabla b_o = (2f)^2$', '$\nabla b_o = (4f)^2$', '$\nabla b_o = (6f)^2$', '$Q_o = 25$ $W m^{-2}$','$Q_o = 100$ $W m^{-2}$','$Q_o = 200$ $W m^{-2}$' , 'Location', 'NorthEastOutside');
% set(l, 'Interpreter' ,'latex', 'FontSize', 20)
% axis equal;
% box on;
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'w', 'Position', [  675         466        1203         508]);

%% ZETA PLOT
figure
hold on
for i=1:nruns;
     if mod(i,4)==1
       set(gca, 'ColorOrderIndex', 1);

     elseif mod(i,4) ==2
                 set(gca, 'ColorOrderIndex', 2);

%          mark = 'o'
     elseif mod(i,4) ==3
                 set(gca, 'ColorOrderIndex', 3);
     else
                 set(gca, 'ColorOrderIndex', 4);

%          mark = 's';
     end
     if i<4
         mark = 'd';
     elseif i>4 & i<9
          mark = 'o';
     else
         mark = 's';
     end
    scatter(abs(imresize(Jbm_zeta(i,:), meanfac)), abs(imresize(Jbs(i,:)+DcoZ.*Jbe(i,:), meanfac)), mark, 'filled', 'MarkerEdgeColor', 'k');
%         scatter(abs(Jbm(i,1:1/meanfac:end)),abs(Jbs(i,1:1/meanfac:end)+Dco.*Jbe(i,1:1/meanfac:end)),mark,  'filled', 'MarkerEdgeColor', 'k');

end
set(gca, 'xscale', 'log', 'yscale', 'log');
set(gca, 'xlim', xl, 'ylim', xl);

xt = linspace(xl(1),xl(end), 100);
plot(xt, xt, 'k');
plot(xt, 2*xt,'--k')
plot(xt, 0.5*xt, '--k')

grid on
