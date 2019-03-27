% calculate Wave Reflection Index using wave intensity analyis

%clear all;
%close all;
%clc;

global P

% Initialisation
addpath(genpath('D:\Documents\Universiteit\Jaar 5\Afstuderen\CircAdaptO212a'));
load('D:\Documents\Universiteit\PhD Jaar 1\Artikel Augmentation Index\CircAdapt\SimulationsTubeArtVen\PNewvMax3kTube12TubeArtVen.mat')
Tube          = 'AoSubclAr';
rho           = 1050; %kg /m3 
Dt            = P.General.Dt;
t             = P.t - P.General.tStart;
tStart        = 1050;
tOk           = tStart : (tStart + (P.General.tCycle / Dt)); 
mmHgToPa      = 133.33;

% % calculation [pTrans, c0, Z0] with anti-collaps stiffness of Aortic -
% % Subclavian segment
% Len   = P.Tube.Len   ;                            % repesentative length of blood vessels
% p0    = P.Tube.p0    ;                            % working pressure
% A0    = P.Tube.A0    ;
% rhob  = P.General.rhob;
% A     = max(1e-10,bsxfun(@rdivide,P.Tube.V,Len)); % vessel cross-section
% ANorm = bsxfun(@rdivide,A,A0);                    % cross-section normalized to physiologic volume
% ap    = 0.02;                                     % introduces: Small volume -> negative transmural pressure pTrans
% mk    = 1 + (P.Tube.k/3-2)/(1+ap);                % stiffness exponential
% pp    = (1+ap) * bsxfun(@power, ANorm, mk);
% pm    = -ap./ANorm;                               % anti-collapse pressure
% pTrans= bsxfun(@times,pp + pm, p0);
% c0    = sqrt( bsxfun(@times, bsxfun(@times,pp,mk)-pm, p0) / rhob );

% Parsing blood pressure
AoBP          = Get('Node','p','Ao'); 
AoBP1cycle    = AoBP(tOk);
AoBP          = AoBP(tOk);
VcBP          = Get('Node','p','Vc'); VcBP = VcBP(tOk);
AoArea        = Get('Tube','A',Tube); AoArea = AoArea(tOk);
AoValveQ      = Get('Valve','q','LvAo'); AoValveQ  = AoValveQ(tOk);
%c0            = c0(:,1); c0 = c0(tOk);

% AoQ is sampled differently -> correct for this observation
AoQ           = Get('Tube','q',Tube); 
AoQ1cycle     = AoQ(tOk);
AoQ           = AoQ(tOk);

% Make first order derivative filters
dP = diff(AoBP) / Dt;
dP(2:end+1) = dP;

Ub = AoQ ./ AoArea;
dU = diff(Ub) / Dt;
dU(2:end+1) = dU;

% Calculate wave speed c0 according to the slope of the PU-relation 
[~ , ~, DnPos] = RelDnHeight(AoBP1cycle);
c0            = CalcWaveSpeedDPDU(dP,dU,DnPos,rho);

% Calculate Wave intensity (dI' = dU/dt . dP/dt) Unit ~= J/m^2 !!!
dI       = dU .* dP;
dIplus   = (dI>0) .* dI;
dImin    = (dI<0) .* dI;

% Wave seperation -> Parker et al. 1990, Khir et al. 2001
dPplus   = 0.5 * (dP + rho .* c0 .* dU) ;
dPmin    = 0.5 * (dP - rho .* c0 .* dU) ;
dUplus   = 0.5 * (dU + (dP ./ (rho*c0))) ;
dUmin    = 0.5 * (dU - (dP ./ (rho*c0))) ;

figure;
plot(dIplus,'b'); hold on;
plot(dImin,'r');
plot(dPplus .* dUplus,'b--')
plot(dPmin .* dUmin,'r--')

Pplus    = (Dt * cumsum(dPplus)) + min(AoBP);
Pmin     = (Dt * cumsum(dPmin))  + min(AoBP);
Uplus    = Dt * cumsum(dUplus);
Umin     = Dt * cumsum(dUmin);

% Find S-wave forward compression wave
SwaveStart = find(dIplus > 1e3,1,'first');
[~,Imax]   = max(dIplus);
SwaveStop  = find(dIplus(Imax:end)==0,1,'first') + Imax;

% Find C-wave reflected decompression wave
CwaveStart = find(dImin ~=0,1,'first');
[~,Imax2]   = min(dImin);
CwaveStop  = find(dImin(Imax2:end)==0,1,'first') + Imax2;

SwaveArea  = Dt * sum(dIplus(SwaveStart:SwaveStop));
CwaveArea  = Dt * sum(dImin(CwaveStart:CwaveStop));

% Calculate Wave Reflection Index, Hughes et al. 2013
WRI = abs(CwaveArea/SwaveArea)
%
figure;
plot(dIplus,'b'); hold on;
plot(dImin,'r');
plot(dPplus .* dUplus,'b--')
plot(dPmin .* dUmin,'r--')
legend('dI','dI>0','dI<0','dP_+ * dU_+','dP-_ * dP_-')
plot(SwaveStart,dIplus(SwaveStart),'ko')
plot(SwaveStop,dIplus(SwaveStop),'ko')
plot(CwaveStart,dImin(CwaveStart),'go')
plot(CwaveStop,dImin(CwaveStop),'go')

figure;
subplot(3,1,1);
plot(AoBP/mmHgToPa,'k'); hold on;
ylabel('Aortic BP [mmHg]')
axis([0,length(dU),min(AoBP/mmHgToPa),max(AoBP/mmHgToPa)])
subplot(3,1,2);
plot(dI,'b'); hold on;
plot(zeros(1,length(dU)),'k','LineWidth',2)
ylabel('Wave Intensity [a.u.]')
axis([0,length(dU),min(dI),max(dI)])
subplot(3,1,3)
plot(sign(dU),'k'); hold on;
plot(sign(dP),'k--')
plot(zeros(1,length(dU)),'k','LineWidth',2)
ylabel('sign of dU and dP [-]')
legend('dU','dP')
axis([0,length(dU),-1.5,1.5])