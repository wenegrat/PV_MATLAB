% foldernames(1,:) = 'GS_025_6F';
% foldernames(2,:) = 'GS_200_6F';
% foldernames(3,:) = 'GS_100_3F';
% foldernames(4,:) = 'GS_200_3F';
% foldernames(5,:) = 'GS_025_2F';%
% foldernames(6,:) = 'GS_200_2F';%
% foldernames(7,:) = 'GS_025_4F';%
% foldernames(8,:) = 'GS_200_4F';

foldernames(1,:) = 'GS_025_2F';
foldernames(2,:) = 'GS_025_4F';
foldernames(3,:) = 'GS_025_6F';
foldernames(4,:) = 'GS_100_2F';
foldernames(5,:) = 'GS_100_4F';%
foldernames(6,:) = 'GS_100_6F';%
foldernames(7,:) = 'GS_200_2F';%
foldernames(8,:) = 'GS_200_4F';
foldernames(9,:) = 'GS_200_6F';


nr = 9;
Jfm = NaN(nr, 1000);
Jbm = Jfm; Jft = Jfm;  Jbt = Jfm;
Jea = Jfm; Jga=Jfm;
for i=1:nr;
    filename = ['./',foldernames(i,:),'/', foldernames(i,:),'_OutputsFlat.mat'];
    load(filename);
    filename = ['./',foldernames(i,:),'/', foldernames(i,:),'_OutputsFlat_2.mat'];
    load(filename);
    nt = length(output.dJf);
    Jfm(i,1:nt) = output.dJf;
    Jbm(i,1:nt) = output.dJb;
%     if i<=4 && i~=3
%     Jft(i,1:nt) = output.dJfga + output.dJfea;
%     Jea(i,1:nt) = output.dJfea;
%     Jga(i,1:nt) = output.dJfga;
%     else
%             Jft(i,1:nt) = output.dJfga - output.dJfea;
%                 Jea(i,1:nt) = -output.dJfea;
%                 Jga(i,1:nt) = output.dJfga;
%     end
%     Jbt(i,1:nt) = output.dJbsa+output.dJbea;
    Jea(i,1:nt) = output2.dJfdavg;
    Jbt(i,1:nt) = output.dJbsa + output2.dJbdavg./ce.*0.06;
    
    cs = strsplit(foldernames(i,:), '_');
    legstring(i,:) = ['Q_o: -', cs{2}, ', M^2: (', cs{3},')^2'];
end

%%
[nruns nt] = size(Jfm);
xl = [1e-5 1e-1];
% H = 150; f0 = 1e-4;Q = 100; gradb = (3*f0).^2;
% B = 9.81*2e-4*Q./(1035*3994);
% Vg = H.*gradb./f0;
% wstar = (abs(B)*H).^(1/3);
% Vttw = .061*wstar.*Vg./(f0.*H);
% Jfi = f0*Vttw*gradb*dx*dy*nx*ny;

meanfac = 1/12;

figure
subplot(1,2,1)
hold on
map = colormap(parula(8));
for i=1:nruns;
     if mod(i,3)==1
       set(gca, 'ColorOrderIndex', 1);

     elseif mod(i,3) ==2
                 set(gca, 'ColorOrderIndex', 2);

%          mark = 'o'
     else
                 set(gca, 'ColorOrderIndex', 3);

%          mark = 's';
     end
     if i<3
         mark = 'd';
     elseif i>3 & i<7
          mark = 'o';
     else
         mark = 's';
     end
%     if i==4; set(gca, 'ColorOrderIndex', 5); end %Skip Purple
    scatter(abs(imresize(Jfm(i,:), meanfac)),0.156./ce.* abs(imresize(Jea(i,:), meanfac)),mark,  'filled', 'MarkerEdgeColor', 'k');

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
tvec = reshape(Jea, nruns*nt, 1);
mask = isfinite(mvec+tvec);
cr = corr(abs(mvec(mask)),abs( tvec(mask)));
title(['Frictional PV Flux,     Corr: ', num2str(cr,2)]);
grid on
xlabel('Model'); ylabel('Scaling');
set(gca, 'FontSize', 16);
subplot(1,2,2)
hold on
for i=1:nruns;
%         if i==4; set(gca, 'ColorOrderIndex', 5); end
     if mod(i,3)==1
       set(gca, 'ColorOrderIndex', 1);

     elseif mod(i,3) ==2
                 set(gca, 'ColorOrderIndex', 2);

%          mark = 'o'
     else
                 set(gca, 'ColorOrderIndex', 3);

%          mark = 's';
     end
     if i<3
         mark = 'd';
     elseif i>3 & i<7
          mark = 'o';
     else
         mark = 's';
     end
    scatter(abs(imresize(Jbm(i,:), meanfac)), abs(imresize(Jbt(i,:), meanfac)), mark, 'filled', 'MarkerEdgeColor', 'k');
    
end
set(gca, 'xscale', 'log', 'yscale', 'log');
set(gca, 'xlim', xl, 'ylim', xl);

xt = linspace(xl(1),xl(end), 100);
plot(xt, xt, 'k');
plot(xt, 2*xt,'--k')
plot(xt, 0.5*xt, '--k')

hold off
mvec = reshape(Jbm, nruns*nt,1);
tvec = reshape(Jbt, nruns*nt, 1);
mask = isfinite(mvec+tvec);
cr = corr(mvec(mask), tvec(mask));
title(['Diabatic PV Flux,     Corr: ', num2str(cr,2)]);
grid on
xlabel('Model'); ylabel('Scaling');
%legend(legstring(1,:), legstring(2,:),legstring(3,:),legstring(4,:), legstring(5,:), legstring(6,:),legstring(7,:),legstring(8,:), legstring(9,:), 'location', 'SouthEastOutside');

% axis equal;
% box on;
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'w', 'Position', [  675         466        1203         508]);

%%
mvec = reshape(Jfm, nruns*nt,1);
evec = reshape(abs(Jea), nruns*nt, 1);
gvec = reshape(abs(Jga), nruns*nt, 1);
mask = isfinite(mvec)
regress(mvec(mask), evec(mask)./ce)

