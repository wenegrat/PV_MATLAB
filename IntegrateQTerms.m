%% INTEGRATE Q TERMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate PV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qi = Q;
Qi(~isfinite(Qi)) = 0;
Qi = Qi.*mask;
Qa = squeeze(nansum(nansum(nansum(Qi.*gridvol)))); %Note this assumes Int_v(Int_t(dQ/dt))
Qt = gradient(Qa, ts); %Take the time derivative
Qa = Qa-Qa(1); % Interested in Delta Q
Qa = Qa./vol; % Average PV substance in the volume (just a normalization factor).