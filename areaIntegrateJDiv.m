%% areaIntegrateJDiv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Easier to integrate the divergence terms for irregular volumes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = time*86400;
%Friction
% Frici = FricDiv;%.*mask;
% Frici(~isfinite(Frici)) = 0;
% Frici = cumtrapz(t, Frici, 4);
% Frici = Frici.*mask;
% Fric = squeeze(trapz(trapz(trapz(Frici)))).*gridvol;
% Fric = Fric./vol;
% 
% Advi = AdvDiv;%.*mask;
% Advi(~isfinite(Advi)) = 0;
% Advi = cumtrapz(t, Advi, 4);
% Advi = Advi.*mask;
% Adv = squeeze(nansum(nansum(nansum(Advi)))).*gridvol;
% Adv = Adv./vol;
% 
% Diai = DiaDiv;%.*mask;
% Diai(~isfinite(Diai))=0;
% Diai = cumtrapz(t, Diai, 4);
% Diai = Diai.*mask;
% Dia = squeeze(nansum(nansum(nansum(Diai)))).*gridvol;
% Dia = Dia./vol;


Frici = FricDiv;
Frici = squeeze(nansum(nansum(nansum(Frici.*mask)))).*gridvol;
Fric = cumtrapz(Frici)*ts;
Fric = Fric./vol;

Diai = DiaDiv;
Diai = squeeze(nansum(nansum(nansum(Diai.*mask)))).*gridvol;
Dia = cumtrapz(Diai)*ts;
Dia = Dia./vol;

Advi = AdvDiv;
Advi = squeeze(nansum(nansum(nansum(Advi.*mask)))).*gridvol;
Adv = cumtrapz(Advi)*ts;
Adv = Adv./vol;
