cd /scratch/jacob13/
% FIRST LOAD ALL OUTPUT FLAT FILES


clear foldernames
% foldernames(1,:) = 'GS_025_2F';
foldernames(1,:) = 'GS_02515F_2D';
foldernames(2,:) = 'GS_02530F_2D';
% foldernames(5,:) = 'GS_100_1F';%
% foldernames(3,:) = 'GS_100_2F';%
% foldernames(4,:) = 'GS_100_4F';%
% foldernames(5,:) = 'GS_100_6F';%
% % foldernames(9,:) = 'GS_200_1F';%
% foldernames(6,:) = 'GS_200_2F';%
% foldernames(7,:) = 'GS_200_4F';
foldernames(3,:) = 'GS_200_6F_2D';

c1 = 2; c2 = 2;

% % foldernames(1,:) = 'GS_025_2F';
% % foldernames(1,:) = 'GS_025_4F';
% % foldernames(2,:) = 'GS_025_6F';
% % foldernames(5,:) = 'GS_100_1F';%
% % foldernames(4,:) = 'GS_100_2F';%
% foldernames(1,:) = 'GS_100_4F';%
% foldernames(2,:) = 'GS_100_6F';%
% % foldernames(9,:) = 'GS_200_1F';%
% % foldernames(7,:) = 'GS_200_2F';%
% foldernames(3,:) = 'GS_200_4F';
% foldernames(4,:) = 'GS_200_6F';
% c1 = 0; c2 = 2;

% foldernames(1,:) = 'GS_025_1F';
% foldernames(2,:) = 'GS_025_2F';
% foldernames(3,:) = 'GS_025_4F';
% foldernames(4,:) = 'GS_025_6F';
% foldernames(5,:) = 'GS_100_1F';%
% foldernames(6,:) = 'GS_100_2F';%
% foldernames(7,:) = 'GS_100_4F';%
% foldernames(8,:) = 'GS_100_6F';%
% foldernames(9,:) = 'GS_200_1F';%
% foldernames(10,:) = 'GS_200_2F';%
% foldernames(11,:) = 'GS_200_4F';
% foldernames(12,:) = 'GS_200_6F';
% 
% c1 = 4; c2 = 8;
[nr, ~] = size(foldernames);

jbsc =1.2; %Empirical coefficient from 1D case.
f0 = 1e-4;
% clear foldernames
% foldernames(1,:) = 'GS_100_4F_LEK';
% nr = 1;

Jfm = NaN(nr, 1000);
Jbm = Jfm; Jft = Jfm;  Jbt = Jfm;
Jea = Jfm; Jga=Jfm;
Jbs = Jfm;
Jbe = Jfm;
TDiff = Jfm;
Etot = TDiff;
CoeffDiff = NaN(nr, 1);
TDiffa = CoeffDiff;
Jbea = CoeffDiff;
Ea = CoeffDiff;
H = Jfm;
De = Jfm;
counter = 0;
meanfac = 1./1;
for i=1:nr;
    try
    filename = ['./',foldernames(i,:),'/', foldernames(i,:),'_OutputsFlat.mat'];
    load(filename);
    

        h = ncread(['./',foldernames(i,:),'/etan.nc'], 'KPPhbl');
   
        nt = length(output.dJf);
        Jfm(i,1:nt) = output.dJf;
        Jbm(i,1:nt) = output.dJb;
        H(i,1:nt)  = squeeze(nanmean(nanmean(h(:,:,:,2:end))));
        FricCoeff = regress(output.dJf, output.dJfea);
        DiaCoeff = regress(output.dJb - jbsc.*output.dJbsa, output.dJbea);
        
        CoeffDiff(i) = (FricCoeff+ DiaCoeff)./(FricCoeff);
%     Jbs(i,1:nt) = output.dJbsa*jbsc;
    avgerr = imresize((output.dJf+ output.dJb - output.dJbsa*jbsc)./(output.dJf), meanfac);
    
    Jbe(i,1:length(avgerr)) = avgerr;
    Jbe(i,1:12) = NaN; % xxx-ad hoc hack
%     Jbe(i,length(avgerr)-48:end) = NaN;
    Jbea(i) = nanmean(Jbe(i,1:nt));
    
    if i<c1+1
        q = 025;
    elseif i>c1 && i<c2+1
        q = 100;
        
    else
       q=200;
    end
     b0 = 9.81*2e-4*squeeze(q)./(1035*3994);
     
     de = squeeze(nanmean(nanmean(sqrt(.1*h.*(b0.*h).^(1/3)./f0))));
     dea = imresize(de, meanfac);     
    
     ha = imresize(squeeze(nanmean(nanmean(h))), meanfac);
     E = squeeze(nanmean(nanmean((.1*h.*(b0.*h).^(1/3)./(h.^2.*f0)))));
     Etot(i,1:length(ha)) = (dea./ha).^2;
%     Etot(i,1:nt) = E(2:end);    
%      TDiff(i,1:nt) = -1+(E(2:end).^(-1/2));
%      TDiffa(i) = nanmean(TDiff(i,1:nt));
         Ea(i) = nanmean((dea./ha).^2);
    H(i,1:length(ha)) = ha;
    De(i,1:length(ha)) = dea;
    catch
        disp(filename)
        counter = counter+1;
    end
    cs = strsplit(foldernames(i,:), '_');
%     legstring(i,:) = ['Q_o: -', cs{2}, ', M^2: (', cs{3},')^2'];
end
disp(['Errors on: ', num2str(counter), ' files.']);
%% FIND BEST FIT COEFFICIENT
[nruns nt] = size(Jfm);
mvec = reshape(Jbe.', nruns*nt,1);
tvec = reshape(TDiff.', nruns*nt, 1);
evec = reshape(Etot.', nruns*nt, 1);
hvec = reshape(H.', nruns*nt, 1);
devec = reshape(De.', nruns*nt, 1);

meanfac = 1./(1);
mveca = imresize(mvec, meanfac);
eveca = imresize(evec, meanfac);
% % subplot(2,1,1)
% scatter(mvec, tvec);
% onetoone
% grid on
% subplot(2,1,2)
% scatter(CoeffDiff, 1./nanmean(H,2));
% % set(gca, 'xlim', [0 1], 'ylim', [0 1])
% grid on 
% onetoone
%%
bindelta = 0.25;
bin = bindelta./2:bindelta:10;
evecb = 30./hvec;sqrt(evec);
evecs = linspace(0,10, 1000);
EV_bin = NaN(size(bin));
e_bin = EV_bin; err_bin = EV_bin;
pts_bin = EV_bin;
err_binn = EV_bin;
for i = 1:length(bin)
    ind = find(evecb>=(bin(i)-bindelta./2) & evecb<(bin(i)+bindelta./2)); % 0.5m/s bins
    pts_bin(i) = length(ind);  % bin nb of points
    if isempty(ind) | length(ind)==1
        EV_bin(i) = NaN; e_bin(i) = NaN; err_bin(i) = NaN; 
    else
        
        EV_bin(i) = nanmedian(mvec(ind)); % bin average power        
        e_bin(i) = bin(i); % bin average speed
        r = quantile(mvec(ind), .75);
        err_bin(i) = r-EV_bin(i); %1.253*2*std(mvec(ind))./sqrt(pts_bin(i)); % bin standard error on the median
        r = quantile(mvec(ind), .25);
        err_binn(i) =  EV_bin(i)-r;
%         
%         bs = bootstrp(1000, @nanmedian, mvec(ind));
%         EV_bin(i) = median(bs); % bin average power        
%         e_bin(i) = bin(i); % bin average speed
%         r = quantile(bs, .025);
%         err_bin(i) = r-EV_bin(i); %1.253*2*std(mvec(ind))./sqrt(pts_bin(i)); % bin standard error on the median
%         r = quantile(bs, .975);
%         err_binn(i) =  EV_bin(i)-r;
    end
end

bfc = regress(mvec, evecb);

ind5 = find(pts_bin>=10);
figure
% plot(sqrt(eveca), mveca, '.', 'MarkerEdgeColor', [1 1 1].*.8)
plot(evecb, mveca, '.', 'MarkerEdgeColor', [1 1 1].*.8)

% plot(bin(ind5), EV_bin(ind5))

hold on
errorbar(e_bin(ind5),EV_bin(ind5),err_binn(ind5), err_bin(ind5), 'LineWidth', 1.5);
% errorbar(e_bin,EV_bin,err_bin, '--');

plot(evecs, evecs, 'k', 'LineWidth', 2);
plot(nanmedian(evecb), (0.2-0.15)./0.2, 'x', 'MarkerSize', 6, 'LineWidth', 2, 'MarkerEdgeColor', 'r'); 
% plot(evecs, evecs.*bfc, '--')
hold off
grid on
xlabel('$\delta_e/\delta_t$', 'FontSize', 16);
ylabel('$ \frac{J_F + J_D - J_{D_{SURF}}}{J_F}$', 'FontSize', 23);
% set(gca, 'XTick', 0:.1:1, 'YTick', 0:.25:1.5);
% set(gca, 'xlim', [0 1], 'ylim', [-.75 2]);
set(gcf, 'Color', 'w' , 'Position', [   675   537   455   437]);
% title(num2str(corr(mvec, evecb, 'rows', 'pairwise')))

%%
figure
mask = isfinite(evec+mvec);
evecb = sqrt(evec) +0;
% binned_plot(evecb(mask), mvec(mask),20)
scatter((evecb), ((evecb).^(-1/4).*15./hvec));
[evecs, isort] = sort(evecb);
evecs = linspace(0,1, 100);
% scatter(evecb, mvec, 'x')
hold on
% for i=1:nruns
%     scatter(Etot(i,:), Jbe(i,:))
% end
% scatter((evecb), ((evecb).^(1/4).*75./hvec));
% scatter(evecb, devec./hvec, 'x')
% plot((evecb), (30./hvec), 'o');
% scatter(evecb, evecb.^(1/2), 'd');
plot(evecs, evecs.^(1/2), 'k', 'LineWidth', 2);
% plot((evecb), ((evecb).^(-1/4).*1./(500*evec)), 'o');
% binned_plot(evecb(mask), mvec(mask),20)
binned_plot(evecb(mask), mvec(mask),20)

hold off
% hold on; plot(nanmedian(evecb), nanmedian(mvec), 'x', 'MarkerSize', 10, 'LineWidth', 4, 'MarkerEdgeColor', 'r'); hold off
hold on; plot(nanmedian(evecb), (0.2-0.15)./0.2, 'x', 'MarkerSize', 10, 'LineWidth', 4, 'MarkerEdgeColor', 'r'); hold off

% set(gca, 'ylim', [-0.5 1.5])
grid on

% %%
% runi = 8;
% subplot(3,1,1)
% scatter( sqrt(Etot(runi,:)), Jbe(runi,:));
% onetoone
% subplot(3,1,2)
% plot(Jbe(runi,:));
% hold on
% plot(sqrt(Etot(runi,:)));
% hold off
% subplot(3,1,3)
% plot(H(runi,:));
% hold on
% plot(De(runi,:))
% hold off
%%
scatter(20./(1*hvec), mvec)
onetoone
grid on
corr(mvec, 20./hvec, 'rows', 'pairwise')
%%
plot((output.dJf + output.dJb-1.2*output.dJbsa)./output.dJf);
hold on
plot(De(end,1:360).'./(4*(squeeze(nanmean(nanmean(h))))), '--');
hold off
