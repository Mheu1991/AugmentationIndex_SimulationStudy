
function [WRI,SwaveArea,CwaveArea] = WaveReflectionIndexFunc(PlotYoN,tStart)

    % calculate Wave Reflection Index using wave intensity analyis

    global P

    % Initialisation
    Tube          = 'AoSubclAr';
    rho           = 1050; %kg /m3 
    Dt            = P.General.Dt;
    t             = P.t - P.General.tStart;
    %tStart        = 1050;
    tOk           = tStart : (tStart + (P.General.tCycle / Dt)); 
    mmHgToPa      = 133.33;
   
    % Parsing blood pressure
    AoBP          = Get('Node','p','Ao'); 
    AoBP1cycle    = AoBP(tOk);
    AoBP          = AoBP(tOk);
    VcBP          = Get('Node','p','Vc'); VcBP = VcBP(tOk);
    AoArea        = Get('Tube','A',Tube); AoArea = AoArea(tOk);
    AoValveQ      = Get('Valve','q','LvAo'); AoValveQ  = AoValveQ(tOk);
    
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
    [~,MaxI ]     = max(AoQ);
    [c0,~,~,~]    = CalcWaveSpeedDPDU(AoBP,Ub,MaxI,rho);
    
    % Calculate Wave intensity (dI' = dU/dt . dP/dt) Unit ~= J/m^2 !!!
    dI = dU .* dP;
    dIplus   = (dI>0) .* dI;
    dImin    = (dI<0) .* dI;

    % Wave seperation -> Parker et al. 1990, Khir et al. 2001
    dPplus   = 0.5 * (dP + rho .* c0 .* dU) ;
    dPmin    = 0.5 * (dP - rho .* c0 .* dU) ;
    dUplus   = 0.5 * (dU + (dP ./ (rho*c0))) ;
    dUmin    = 0.5 * (dU - (dP ./ (rho*c0))) ;

    Pplus    = (Dt * cumsum(dPplus)) + min(AoBP);
    Pmin     = (Dt * cumsum(dPmin))  + min(AoBP);
    Uplus    = Dt * cumsum(dUplus);
    Umin     = Dt * cumsum(dUmin);

    % Find S-wave forward compression wave
    SwaveStart = find(dIplus > 1e3,1,'first');
    [~,Imax]   = max(dIplus);
    SwaveStop  = find(dIplus(Imax:end)==0,1,'first') + Imax;

    % Find C-wave reflected decompression wave
    CwaveStart = find(dImin<-1000,1,'first');
    [~,Imax2]   = min(dImin);
    CwaveStop  = find(dImin(Imax2:end)==0,1,'first') + Imax2;

    SwaveArea  = Dt * sum(dIplus(SwaveStart:SwaveStop));
    CwaveArea  = Dt * sum(dImin(CwaveStart:CwaveStop));

    % Calculate Wave Reflection Index, Hughes et al. 2013
    WRI = abs(CwaveArea/SwaveArea);
    CwaveArea = abs(CwaveArea);
    %
    
    %% Plotting
    if strcmp(PlotYoN,'true')
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
        axis([0,900,-1e6,2.5e6])
    end
%%
    
end

