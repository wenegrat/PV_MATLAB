%% INTEGRATE Q TERMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate PV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Qi = Q+cumsum(dADV, 4).*ts; % XX-TESTING mass weighted...
Qi = Q;%+cumsum(dADV,4).*ts; 
Qi(~isfinite(Qi)) = 0;
% rho = 1035*( 1-2e-4*(THETA - 16));
% maskt = mask;
Qi = Qi.*mask;
Qa = squeeze(nansum(nansum(nansum(Qi.*gridvol)))); %Note this assumes Int_v(Int_t(dQ/dt))
% Qa = squeeze(trapz(Z, trapz(Y, trapz(X, Qi))));
Qt = gradient(Qa, ts); %Take the time derivative
Qa = Qa-Qa(1); % Interested in Delta Q
Qa = Qa./vol; % Average PV substance in the volume (just a normalization factor).