%% Calculate and plot fibre strain rate 
% M Heusinkveld 21-08-2018
function PeakSR = StrainRateFunc

global P

Efiber = exp((P.Patch.Ef(:,4)+P.Patch.Ef(:,4))/2); % logarithmic to cauchy strain
SRfiber= gradient(Efiber')/P.General.Dt;

PeakSR = abs(min(SRfiber));



