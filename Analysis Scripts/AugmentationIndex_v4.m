function [AIx,TTP,AIxRev,iInflection,AugmentationPressure] = AugmentationIndex_v4(Pressure,iS)  

% This function provides a tool to determine the Augmentation Index of a
% Aortic BP waveform
% Window to calculate 2nd-order derivative is pAO (pAO > MAP & pAO <= Psys)
% Edited Maarten Heusinkveld 22-10-2014

global P

t    = P.t ;
t    = t - t(1);
dt   = P.General.Dt;

dk = P.Tube.k(1) - 8;


if nargin <2
    iS  = 1050;   
end

% Determine the 2nd derivative of the aortic bp waveform
% df/dt = [f(x+dt)-f(x)]/dt

d1S = diff(Pressure)/dt;
D1S(2:length(Pressure)) = d1S; D1S(1) = D1S(2);

d2S = diff(D1S)/dt; 
D2S(2:length(D1S))  = d2S; D2S(1) = D2S(2);

% Select a single/ half beat, assuming steady state
tOK  = iS:(iS + P.General.tCycle/dt);
tOKh = iS-20:(iS + P.General.tCycle/(2*dt));
stOK = length(tOK); 

D2Temp   = D2S(tOK);
AObpTemp = Pressure(tOK);
MeanBP   = mean(Pressure);
MaxBP    = max(AObpTemp);
MinBP    = min(AObpTemp);

% Augmentation is assumed to be above MAP and under Systolic Pressure: store cutoff index
%IdxStart        = find(AObpTemp>0.5*MeanBP,1,'first');

if dk <= -1
    IdxStart = 100;
else
    IdxStart = 70;
end

IdxStop         = find(AObpTemp == MaxBP,1,'first');

%%
IdxMin          = find( D2S(tOKh) == max(D2S(tOKh)));

%%
D2Augmented     = D2Temp(IdxStart:IdxStop);
%[~,I]           = max(D2Augmented);
%plot(D2Augmented)
% Use sign change to calculate First zero crossing index corresponds with f
iInflection = find(D2Augmented==0);
if isempty(iInflection)
  % At the moment of crossing, the sign will change:
  % find zero crossing before maximum 2nd-derivative value is achieved
  s = diff(sign( D2Augmented ) );
  iInflection = find(s,1,'first') + IdxStart;
  if isempty(iInflection)
      AIx = NaN;
      AIxRev = NaN;
      TTP = NaN;
      AugmentationPressure = NaN;
      return
  end
end

% Calculate corresponding augmentation pressure (AP) and calculate
% Augmentation Index ,AIx = (SBP-AP)/PP (*100)%
AugmentationPressure = AObpTemp( iInflection );
PeakSystP            = max(AObpTemp);
DiastP               = min(AObpTemp);
PP                   = PeakSystP - DiastP;
AIx                  = (PeakSystP - AugmentationPressure)/ PP;
AIxRev               = (PeakSystP - AugmentationPressure) / (AugmentationPressure - DiastP);
TTP                  = t(iInflection) - t(IdxMin) ;

% if nargin > 3
%     figure;
%     subplot(2,1,1); plot(AObpTemp); hold on
%     plot(iInflection,AObpTemp( iInflection ),'ro')
%     axis([0,350,40,180])
%     subplot(2,1,2); plot(D2Temp); hold on;
%     plot(tOK-tOK(1),zeros(1,stOK),'r');
%     plot(tOK-tOK(1),mean(AObpTemp*ones(1,stOK)),'r--');
%     plot(iInflection,D2Temp( iInflection ),'ro')
%     axis([0,350,-1e4,1e4])
%     disp(['AugmentationIndex = ' num2str(100*AIx) ' %'])
% end

end
