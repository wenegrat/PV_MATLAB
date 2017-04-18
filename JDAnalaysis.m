[~, ~, ~, bt] = gradient(b, ts);

% bt = TtoB.*GetVar(statefile, diagfile, {'TOTTTEND', '(1)./86400'}, slice);


Jb = outputFull.JBz(:,:,2,st:(st+nt-1));
Jba = squeeze(nanmean(nanmean(Jb)));


bta = f0.*squeeze(nanmean(nanmean(bt(:,:,2,:))));

bzet = squeeze(nanmean(nanmean(bt(:,:,2,:).*ZETA(:,:,2,:))));

bad = U.*bx + V.*by;
bad = bad(:,:,2,:);
bad = squeeze(nanmean(nanmean( bad(:,:,1,:).*(f0+ZETA(:,:,2,:)))));

wterms = squeeze(nanmean(nanmean( (f0+ZETA(:,:,2,:)).*W(:,:,2,:).*bz(:,:,2,:))));


wpbp = Drv(metric, W.*b, 'z');
%%
wpbar = squeeze(nanmean(nanmean(wpbp(:,:,2,:).*(f0))));
wpbarz = squeeze(nanmean(nanmean(wpbp(:,:,2,:).*(ZETA(:,:,2,:)))));

%%
figure
subplot(2,1,1)
plot(Jba)
hold on
plot(-bta);
plot(-bzet);
plot(-bad);
plot(-wterms);
plot(-bta - bzet-bad-wterms,'--');
plot(-wpbar,'r')
% plot(-bta -bzet -wpbar,'g')
hold off
legend('$J_D$', '$f b_t$', '$\zeta b_t$','$U_h \nabla b$', '$WB_z$' )
subplot(2,1,2)
scatter(Jba, -bta-bzet-bad-wterms)
onetoone
