% One dimensional budget

xi = 40; 
yi = 30;
QT = NaN(size(squeeze(Q(xi,yi,:,:))));
QT(:,2:end-1) =squeeze(( Q(xi,yi,:,3:end) - Q(xi,yi,:,1:end-2))./(2*ts));
% q = squeeze(QT(:,ti));
dx = 500;
dy = 500;%Approx.
dz1 = repmat(zm, [1 nt]);

uq = Uf.*squeeze(Q(:,:,:,:));
vq = Vf.*squeeze(Q(:,:,:,:));
wq = Wf.*squeeze(Q(:,:,:,:));
duqdx = squeeze((uq(xi+1,yi,:,:)-uq(xi-1,yi,:,:))./(2*dx));
dvqdy = squeeze((vq(xi,yi+1,:,:)-vq(xi,yi-1,:,:))./(2*dy));
dwqdz = NaN(size(dvqdy));
dwqdz(2:end-1,:) = squeeze(squeeze((wq(xi,yi,3:end,:) - wq(xi,yi,1:end-2,:)))./((dz1(3:end,:) - dz1(1:end-2,:))));

ADV = (duqdx+dvqdy+dwqdz);
%%
ti = 20;
q = squeeze(QT(:,ti));
qtot = q+ADV(:,ti);
% qtots = smooth(qtot,3);
figure
plot(q, zm)

hold on
plot(-ADV(:,ti), zm)
plot(qtot,zm, '--')
% plot(-dwqdz(:,ti), zm)
% plot(-(duqdx(:,ti) + dvqdy(:,ti)), zm)
% plot(-squeeze(JFT(xi,yi,ti)), 0, 'x')
% plot(-squeeze(JBT(xi,yi,ti)), 0, 'x')
% plot(-squeeze(JFW(xi,yi,ti)), 0, 'x')
% plot(-squeeze(JBT(xi,yi,ti)+JFT(xi,yi,ti)+JFW(xi,yi,ti)), 0, 'x')
% plot(cumtrapz(squeeze(z(xi,yi,3:end-2)),qtots(3:end-2)), zm(3:end-2))
% xt = get(gca, 'xtick');
% plot(xt, squeeze(-hkpp(xi,yi,ti)).*ones(size(xt)));
hold off
set(gca, 'ylim',[-1000 0])
grid on
%%
% scatter(qtot, -dwqdz(:,ti))
%%
zi = 28;

figure
plot(QT(zi,:));
hold on
plot(-ADV(zi,:));
hold off

%%
figure
scaletot = JFT + JBT +JFW;
% scaletot = JFT+JBT;
scaletot = squeeze(scaletot(xi,yi,:));
Qtot = QT+ADV;
Qnc = NaN(nt,1);
for i=1:nt
    dz2 = squeeze(dz(xi,yi,:,i));
    mask = isfinite(Qtot(:,i));
%     mask = mask & zm >= -squeeze(hkpp(xi,yi,i));
Qnc(i) = trapz(dz2(mask), Qtot(mask,i));
end
plot(Qnc)
hold on
plot(-scaletot)
hold off

%%
figure
plot(cumtrapz(Qnc));
hold on
plot(cumtrapz(-scaletot));
    hold off

mask = isfinite(Qnc+scaletot);
corr(cumtrapz(Qnc), cumtrapz(-scaletot), 'rows', 'pairwise')