function [Dn , f , absf] = RelDnHeight(AObp)
%Determine the absolute position of the dicrotic notch at the aortic bp
%curve by identifying the maximum of the second-derivative

global P

dt   = P.General.Dt;

if nargin < 1
    
    AObp = Get('Node','p','Ao')/133;
    
    % Determine the 2nd derivative of the aortic bp waveform
    % df/dt = [f(x+dt)-f(x)]/dt

    d1S = diff(AObp)/dt;
    D1S(2:length(AObp)) = d1S; D1S(1) = D1S(2);

    d2S = diff(D1S)/dt; 
    D2S(2:length(D1S))  = d2S; D2S(1) = D2S(2);

    % Select a single beat, assuming steady state
    tOK = P.General.tCycle/dt;

    D2Temp      = D2S(1:tOK);
    AObpTemp    = AObp(1:tOK);
    
else
    d1S = diff(AObp)/dt;
    D1S(2:length(AObp)) = d1S; D1S(1) = D1S(2);

    d2S = diff(D1S)/dt; 
    D2S(2:length(D1S))  = d2S; D2S(1) = D2S(2);

    
    D2Temp      = D2S;
    AObpTemp    = AObp;
end

% Dn occurs after systolic BP is reached
IdxSys      = find(AObpTemp == max(AObpTemp));
AObpEStoDIA = AObpTemp( IdxSys : end-100);
D2EStoDIA   = sgolayfilt(D2Temp(IdxSys : end-100),4,21);


f    = find(D2EStoDIA==max(D2EStoDIA));
absf = IdxSys + f;
Dn   = (AObpEStoDIA(f) - min(AObpTemp))/ (max(AObpTemp)-min(AObpTemp));

end

