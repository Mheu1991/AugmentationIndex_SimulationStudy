% Calculate wave intensity from aortic pressure and flow across the aortic
% valve
clear all;
close all;
clc;

addpath('D:\Documents\Universiteit\PhD Jaar 1\Artikel Augmentation Index\CircAdapt\');

load('D:\Documents\Universiteit\PhD Jaar 1\Artikel Augmentation Index\CircAdapt\SimulationsTubeArtVen\PNewvMax7kTube12TubeArtVen.mat')

% type of analysis wave seperation / reservoir wave
%Analysis = 'wavesep';
Analysis = 'reswave';
%
Tube          = 'AoSubclAr';
%Tube          = 'AoBrAr';
%
rho           = 1050; %kg /m3 
Dt            = P.General.Dt;
t             = P.t - P.General.tStart;

tStart        = 1050;
tOk           = tStart : (tStart + (P.General.tCycle / Dt)); 

%AoBP          = Get('Node','p','Ao'); 
AoBP          = Get('Tube','p',Tube); 
AoBP1cycle    = AoBP(tOk);
AoBP          = AoBP(tOk);
VcBP          = Get('Node','p','Vc'); VcBP = VcBP(tOk);
AoArea        = Get('Tube','A',Tube); AoArea = AoArea(tOk);
AoValveQ      = Get('Valve','q','LvAo'); AoValveQ  = AoValveQ(tOk);

%+++%
% c0            = Get('Tube','c0',Tube);
% c0            = c0(tOk);
% c0            = rms(c0); 
%+++%

t1cycle       = t(tOk); % for Windkessel calculation
t             = t(tOk);
mmHgToPa      = 133.33;

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
[~ , ~, DnPos]      = RelDnHeight(AoBP1cycle);
[~,MaxI ]           = max(AoQ);
[c0,P_s,U_s,ftk]    = CalcWaveSpeedDPDU(AoBP,Ub,MaxI,rho);
c0
%c0 = 20;
% Calculate Wave intensity (dI' = dU/dt . dP/dt) Unit ~= J/m^2 !!!
dI = dU .* dP;
%

switch Analysis
    
    case 'wavesep'
        % Calculate wave intensity and pressure contribution of forward waves and backwards
        % waves
        dIplus   = (dI>0) .* dI;
        dImin    = (dI<0) .* dI;
                
        dPplus   = 0.5 * (dP + rho .* c0 .* dU) ;
        dPmin    = 0.5 * (dP - rho .* c0 .* dU) ;
        dUplus   = 0.5 * (dU + (dP ./ (rho*c0))) ;
        dUmin    = 0.5 * (dU - (dP ./ (rho*c0))) ;
                
        Pplus    = (Dt * cumsum(dPplus)) + min(AoBP);
        Pmin     = (Dt * cumsum(dPmin))  + min(AoBP);
        Uplus    = Dt * cumsum(dUplus);
        Umin     = Dt * cumsum(dUmin);
        
        pb = max(abs(Pplus));
        pf = max(abs(Pmin));
        pbpf = pb / pf
        
        % Compare (normalized) forward wave pressure with aortic flow
        % waveforms
        
        PplusN =  (Pplus - min(Pplus)) ./ range(Pplus);
        AoQN =  (AoQ - min(AoQ)) ./ range(AoQ);       
        
        
        % AUC of dIplus and dImin
        AUCdIplus            = trapz(t,dIplus);
        AUCdImin             = trapz(t,dImin);
        RatioPlusMin         = abs(AUCdImin/AUCdIplus);
        
        %% Plotting      

        P0 = min(AoBP)/133;

        figure;
        plot(AoBP/133,'c--','LineWidth',2); hold on;
        plot( ((Pmin-Pmin(1))/133)+((Pplus-Pplus(1))/133) + P0 ,'k','LineWidth',1); 
        plot(((Pmin-Pmin(1))/133),'b','LineWidth',2); 
        plot(((Pplus-Pplus(1))/133),'r','LineWidth',2)
        plot(zeros(1,length(Pplus)),'k')
        axis([0 900 -30 180])
        ylabel('pressure [mmHg]')
        xlabel('time [ms]')
        legend('P_s','P_f + p_b ','P_b','P_f')
        
        figure;
        plot(Ub,'c--','LineWidth',2); hold on;
        plot( Umin+Uplus ,'k','LineWidth',1); 
        plot(Umin,'b','LineWidth',2); 
        plot(Uplus,'r','LineWidth',2)
        plot(zeros(1,length(Pplus)),'k')
        legend('U_s','U_f + U_b ','U_b','U_f')
        
        figure;
        plot(dU_dia,dP_sys)
        xlabel('U')
        ylabel('P')

        
    case 'reswave'
                      
        figure;
        [CFitted,RpFitted,ResP,P0F,resPerc] = FindReservoirPressure(AoBP,AoBP1cycle,AoQ1cycle,t1cycle,Dt);
        Pres                                = WKfunc(CFitted,RpFitted,t1cycle,Dt,P0F,AoQ * 1e6) ; % Pa
        Pwave                               = AoBP' - Pres;         
        
        % Make first order derivative filters
        dP          = diff(Pwave) / Dt;
        dP(2:end+1) = dP;
        dU          = diff(Ub) / Dt;
        dU(2:end+1) = dU;

        % Calculate Wave intensity (dI' = dU/dt . dP/dt) Unit ~= J/m^2 !!!
        dI          = dU .* dP';
        
        % Calculate wave intensity and pressure contribution of forward waves and backwards
        % waves
        dIplus = (dI>0) .* dI;
        dImin  = (dI<0) .* dI;

        dPplus = 0.5 * (dP' + rho .* c0 .* dU) ;
        dPmin  = 0.5 * (dP' - rho .* c0 .* dU) ;
        dUplus = 0.5 * (dU + (dP' ./ (rho*c0))) ;
        dUmin  = 0.5 * (dU - (dP' ./ (rho*c0))) ;

        Pplus  = (Dt * cumsum(dPplus)) + min(Pwave) ;
        Pmin   = (Dt * cumsum(dPmin))  + min(Pwave) ;
        Uplus  =  Dt * cumsum(dUplus);
        Umin   =  Dt * cumsum(dUmin);
        
        pb = max(abs(Pmin));
        pf = max(abs(Pplus));
        pbpf = pb/ pf
        
        % calculate summed signal for comparison purposes with simulated
        % signal
        P0 = min(AoBP)/mmHgToPa;                        %diastolic pressure
        pbpfpr = ((Pmin-Pmin(1))/mmHgToPa)+((Pplus-Pplus(1))/mmHgToPa) + (Pres'-Pres(1))/mmHgToPa + P0;
        
        % Compare (normalized) forward wave pressure with aortic velocity
        % waveforms
        
        PplusN = (Pplus - min(Pplus)) ./ range(Pplus);
        AoQN   = (AoQ - min(AoQ)) ./ range(AoQ);       
        
        %% Plotting   

        figure;
        plot(AoBP/mmHgToPa,'c--','LineWidth',2); hold on;
        plot( pbpfpr ,'k','LineWidth',1); 
        plot(Pwave/mmHgToPa,'k--','LineWidth',2)
        plot(((Pmin-Pmin(1))/mmHgToPa),'b','LineWidth',2); 
        plot(((Pplus-Pplus(1))/mmHgToPa),'r','LineWidth',2)
        plot((Pres - Pres(1))/mmHgToPa,'g')
        plot(zeros(1,length(Pres)),'k')
        legend('p_s','P_f + p_b + p_r','p_{wave}','p_b','p_f','p_{res}')
        ylabel('pressure [mmHg]')
        xlabel('time [ms]')
        axis([0 900 -30 180])
        
        figure;
        plot(Ub,'c--','LineWidth',2); hold on;
        plot( Umin+Uplus ,'k','LineWidth',1); 
        plot(Umin,'b','LineWidth',2); 
        plot(Uplus,'r','LineWidth',2)
        plot(zeros(1,length(Pres)),'k')
        legend('U_s','U_f + U_b ','U_b','U_f')
        
        figure;
        plot(Ub,AoBP/mmHgToPa,'b'); hold on;
        plot(U_s,P_s/mmHgToPa,'r')
        FitL = ftk.a1 + (ftk.a2 .* U_s);
        plot(U_s,FitL/mmHgToPa,'k','LineWidth',2)
        xlabel('U [m/s]')
        ylabel('P [mmHg]')
        legend('Total PU loop','PU loop during systole','Fit of the linear portion')
        

end