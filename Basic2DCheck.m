load('GS_025SMF_2D_OutputsFull.mat')
load('GS_025SMF_2D_OutputsFlat.mat')

%%
plot(output.Qa);
hold on
plot(-(output.Jfa+output.Jba),'--')
plot(-cumtrapz(abs(output.dJfea).*0.025).*7200./squeeze(nansum(nansum(nansum(outputFull.gridvol)))))
hold off
%%

%%
plot(smooth(output.Qt, 1))
hold on
% plot(-output.dJf)
% plot(.15.*output.dJfea, '--')
% plot(-output.dJb)
% plot(-.118.*output.dJbea, '--')
% plot(.031.*output.dJfea, '--')

plot(-smooth(output.dJf+output.dJb,1))
plot(-output.dJbsa)
plot(output.dJfea.*0.013, '--')
hold off
grid on

%%
tl = 13:length(output.dJf);
regress(output.dJf(tl), output.dJfea(tl))
% scatter(output.dJf(tl), output.dJfea(tl))
regress(output.dJb(tl)-1.2.*output.dJbsa(tl), output.dJbea(tl))

% corr(output.dJf(tl), output.dJfea(tl))
%%
ustar = (output.H.*9.81*2e-4*squeeze(25)./(1035*3994)).^(1/3);
Coeff = regress(output.dJf, ustar./1e-4.*output.MagB.^2)
scatter(output.dJf, Coeff.*ustar./1e-4.*output.MagB.^2)
corr(output.dJf, Coeff.*ustar./1e-4.*output.MagB.^2)

grid on
onetoone