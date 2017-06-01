%% 1) Go To folder of interest
%% 2) Load OutputsFull
%% 3) Run this stuff
statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc'; extrafile = 'extra.nc';
kppfile = 'kppdiags.nc';

TtoB = 9.81.*2e-4;
f0 = 1e-4;

Q0 = ncread(etanfile, 'TFLUX');
X = ncread(statefile, 'X');
Y = ncread(statefile, 'Y');
Z = ncread(statefile, 'Z');
Zl = ncread(statefile, 'Zl');
T = ncread(diagfile, 'T');

nx = length(X); ny = length(Y); nz = length(Z);
try
dx = X(2)-X(1)
catch
    dx = 500;
end
dy = Y(2)-Y(1)
dz = Z(1)-Z(2) %surface only, XX-should track this through the code to ensure correct.
ts = T(2)-T(1)

dh = diff([Zl; -300]);

tslice = [1 length(T)-1];
nt = tslice(end);
gridvol = permute(repmat( dx.*dy.*abs(dh), [1 nx ny nt]), [2 3 1 4]);

    

slice={0, 0, 0, tslice};
sliceEta={0,0,[1 1],tslice};

%%
zetas = squeeze(GetVar(statefile, extrafile, {'momVort3', '(1)'}, sliceEta))+f0;

JBN = squeeze(outputFull.JBz(:,:,2,:)-outputFull.JBz(:,:,end-1,:));
%%
zetasvec = reshape(zetas./f0, [nx.*ny*nt, 1]);
jbvec = reshape(JBN.*dx.*dy, [nx.*ny*nt, 1]);
edges = -2:0.5:10;
[N,edges,bin] = histcounts(zetasvec, edges);
    
m = bin>0;
AvgJB = accumarray(bin(m),jbvec(m),[],@nansum);
edgesa = 0.5.*(edges(2:end)+edges(1:end-1));
AvgJB = AvgJB(AvgJB~=0);
m = (N>0);
edgesa = edgesa(m);
%% HISTOGRAM FIGURE
figure
% subplot(1,2,1);
bar(edgesa, -AvgJB)
ylabel('$-\int J_D^{NUM} \mathrm{d A}$   $\left[m^3s^{-4}\right]$')
set(gca, 'FontSize', 20, 'XTick', -2:2:10)
xlabel('$\frac{f+\zeta}{f}$', 'FontSize', 24);

set(gca, 'xlim', [-2 10]);

grid on
set(gcf, 'Color', 'w')

%% TRY REGRESSION
jbna = squeeze(nanmean(nanmean(JBN)));
jbnaz = squeeze(nanmean(nanmean(f0./(zetas).*JBN)));

jbtae = output.dJbea./(dx.*dy*nx.*ny);
jbtas = output.dJbsa./(dx.*dy.*nx.*ny);

regress(jbna-1.*1.2*jbtas,abs( jbtae))
scatter(jbna-1.*1.2*jbtas, jbtae)
regress(jbnaz-1.*1.2*jbtas, jbtae)

