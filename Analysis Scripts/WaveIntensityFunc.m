function [pbpf,pbpfpb,CFitted,RpFitted,resPerc,pb,pf,c0] = WaveIntensityFunc(Analysis,PlotYoN,tStart)
    % Calculate wave intensity from aortic pressure and flow across the aortic
    % valve
     
    global P

    % type of analysis wave seperation (see Parker1990) / reservoir wave
    % (see Wang2003) / fourier (see Westerhof 1972)
    
    Tube          = 'AoSubclAr';
    rho           = 1050; %kg /m3 
    mmHgToPa      = 133.33;
    
    Dt            = P.General.Dt;
    sf            = 1/Dt;
    freq          = 1/P.General.tCycle;
    
    t             = P.t - P.General.tStart;
    tOk           = tStart : (tStart + (P.General.tCycle / Dt));
    t1cycle       = t(tOk); % for Windkessel calculation
    t             = t(tOk);
    
    AoBP          = Get('Node','p','Ao'); 
    AoBP1cycle    = AoBP(tOk);    
    AoBP          = AoBP(tOk);
    pd            = min(AoBP);
    VcBP          = Get('Node','p','Vc'); VcBP = VcBP(tOk);
    
    AoArea        = Get('Tube','A',Tube); AoArea = AoArea(tOk);
    AoValveQ      = Get('Valve','q','LvAo'); AoValveQ  = AoValveQ(tOk);
    AoQ           = Get('Tube','q',Tube); 
    AoQ1cycle     = AoQ(tOk);
    AoQ           = AoQ(tOk);

    % Make first order derivative filters
    dP = diff(AoBP) / Dt;
    dP(2:end+1) = dP;
    
    % Calculate the blood flow velocity and the incremental change
    Ub = AoQ ./ AoArea;
    dU = diff(Ub) / Dt;
    dU(2:end+1) = dU;
    
    % Calculate wave speed c0 according to the slope of the PU-relation 
    [~ , ~, DnPos] = RelDnHeight(AoBP1cycle);
    [~,MaxI]       = max(AoQ);
    c0             = CalcWaveSpeedDPDU(AoBP,Ub,MaxI,rho);

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

            Pplus    = (Dt * cumsum(dPplus));% + pd;
            Pmin     = (Dt * cumsum(dPmin));%  + pd;
            Uplus    = Dt * cumsum(dUplus);%   + 0;
            Umin     = Dt * cumsum(dUmin) ;%   + 0;
            
            %+++%
            [pb,imax_pb] = max(abs(Pmin));
            pb           = (pb)/mmHgToPa;
            [pf,imax_pf] = max(abs(Pplus));
            pf           = (pf)/mmHgToPa;
            pbpf         = pb /pf;
            pbpfpb       = pb / (pb+pf);
            
            
            CFitted  = NaN;  %Auxiliary
            RpFitted = NaN;
            resPerc  = NaN;
            %RefCoeff2= max(Get('Tube','pP','AoSubclAr')) / max(Get('Tube','pD','AoSubclAr'));
           
            % Compare (normalized) forward wave pressure with aortic flow
            % waveforms

                        
            if strcmp(PlotYoN,'true')
                
                P0 = min(AoBP)/133;
                figure;
                plot(AoBP/133,'c--','LineWidth',2); hold on;
                plot( ((Pmin-Pmin(1))/mmHgToPa)+((Pplus-Pplus(1))/mmHgToPa) + P0 ,'k','LineWidth',1); 
                plot(((Pmin)/mmHgToPa),'b','LineWidth',2); 
                plot(((Pplus)/mmHgToPa),'r','LineWidth',2)
                plot(zeros(1,length(Pplus)),'k')
                plot(imax_pb,(Pmin(imax_pb))/mmHgToPa,'ko')
                plot(imax_pf,(Pplus(imax_pf))/mmHgToPa,'ko')
                axis([0 900 -30 180])
                ylabel('pressure [mmHg]')
                xlabel('time [ms]')
                legend('P_s','P_f + p_b ','P_b','P_f')

%                 figure;
%                 plot(Ub,'c--','LineWidth',2); hold on;
%                 plot( Umin+Uplus ,'k','LineWidth',1); 
%                 plot(Umin,'b','LineWidth',2); 
%                 plot(Uplus,'r','LineWidth',2)
%                 plot(zeros(1,length(Pplus)),'k')
%                 legend('U_s','U_f + U_b ','U_b','U_f')
% 
%                 figure;
%                 plot(dU_dia,dP_sys)
%                 xlabel('U')
%                 ylabel('P') 
                
            end
          

        case 'reswave'
                       
            [CFitted,RpFitted,~,P0F,resPerc] = FindReservoirPressure(AoBP,AoBP1cycle,AoQ1cycle,t1cycle,Dt);
            Pres  = WKfunc(CFitted,RpFitted,t,Dt,P0F,AoQ * 1e6) ; % Pa
            Pwave = AoBP' - Pres;         

            % Make first order derivative filters
            dP = diff(Pwave) / Dt;
            dP(2:end+1) = dP;

            Ub = AoQ ./ AoArea;
            dU = diff(Ub) / Dt;
            dU(2:end+1) = dU;

            % Calculate Wave intensity (dI' = dU/dt . dP/dt) Unit ~= J/m^2 !!!
            dI = dU .* dP';

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

            [pb,imax_pb]     = max(abs(Pmin));
            [pf,imax_pf]     = max(abs(Pplus));
            pb     = pb/mmHgToPa;
            pf     = pf/mmHgToPa;
            pbpf   = pb/(pf+pd);
            pbpfpb = pb / (pb+pf);
            
            % calculate summed signal for comparison purposes with simulated
            % signal
            P0     = min(AoBP)/mmHgToPa;                        %diastolic pressure
            pbpfpr = ((Pmin-Pmin(1))/mmHgToPa)+((Pplus-Pplus(1))/mmHgToPa) + (Pres'-Pres(1))/mmHgToPa + P0;
            

            %%
            
            if strcmp(PlotYoN,'true')
                figure;
                plot(AoBP/mmHgToPa,'c--','LineWidth',2); hold on;
                plot( pbpfpr ,'k','LineWidth',1); 
                plot(Pwave/mmHgToPa,'k--','LineWidth',2)
                plot(Pmin/mmHgToPa,'b','LineWidth',2); 
                plot(Pplus/mmHgToPa,'r','LineWidth',2)
                plot((Pres - Pres(1))/mmHgToPa,'g')
                plot(zeros(1,length(Pres)),'k')
                plot(imax_pb,Pmin(imax_pb)/mmHgToPa,'ko')
                plot(imax_pf,Pplus(imax_pf)/mmHgToPa,'ko')
                legend('p_s','P_f + p_b + p_r','p_{wave}','p_b','p_f','p_{res}')
                ylabel('pressure [mmHg]')
                xlabel('time [ms]')
                axis([0 900 -30 180])

%                 figure;
%                 plot(Ub,'c--','LineWidth',2); hold on;
%                 plot( Umin+Uplus ,'k','LineWidth',1); 
%                 plot(Umin,'b','LineWidth',2); 
%                 plot(Uplus,'r','LineWidth',2)
%                 plot(zeros(1,length(Pres)),'k')
%                 legend('U_s','U_f + U_b ','U_b','U_f')
% 
%                 figure;
%                 plot(Ub,AoBP/mmHgToPa,'b'); hold on;
%                 plot(U_s,P_s/mmHgToPa,'r')
%                 FitL = ftk.a1 + (ftk.a2 .* U_s);
%                 plot(U_s,FitL/mmHgToPa,'k','LineWidth',2)
%                 xlabel('U [m/s]')
%                 ylabel('P [mmHg]')
%                 legend('Total PU loop','PU loop during systole','Fit of the linear portion') 
                
            end           
        case 'fourier'
        
                noHarmonics= 15;
                ns         = 1e3;
                Pcoef      = fourier(AoBP,sf,noHarmonics,freq);
                Qcoef      = fourier(AoQ,sf,noHarmonics,freq);
                Zcoef      = Pcoef ./ Qcoef;
                ZcoefAngle = angle(Zcoef);

                %
                % Seperate waveform in forward and backward component according to
                % Westerhof et al. 1972 )

                % Zc: characteristic impedance is taken as the minimal of Zin in the 0-10
                % Hz interval

                Zc             = mean(Zcoef(end-5:end));

                % calculate the moduli of pf and pb
                for n = 1: noHarmonics+1

                    Z          = Zcoef(n);        

                    Rc(n)      = (Z - Zc) / (Z + Zc);    

                    PFcoef(n)      = Pcoef(n) /(1+Rc(n));
                    PBcoef(n)      = (Pcoef(n) * Rc(n)) / (1+Rc(n));

                    QFcoef(n)      = Qcoef(n) /(1-Rc(n));
                    QBcoef(n)      = -Rc(n) * (Qcoef(n) / (1-Rc(n))); 

                end

                % inverse fourier transformation to time-domain signal

                Pplus = ifourier(PFcoef, ns, freq);
                Pmin  = ifourier(PBcoef, ns, freq);
                Qplus = ifourier(QFcoef, ns, freq);
                Qmin  = ifourier(QBcoef, ns, freq);

                % calculate relection index amplitude(pb)/amplitude(pf)

                pb = max(Pmin)/ mmHgToPa;
                pf = max(Pplus)/ mmHgToPa;

                pbpf     = pb/pf;
                pbpfpb   = pb/ (pf+pb);
                CFitted  = NaN;  %Auxiliary
                RpFitted = NaN;
                ErrFitted= NaN;
                resPerc  = NaN;
                c0       = NaN;
    end
end