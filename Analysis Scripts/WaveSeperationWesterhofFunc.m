
function [pb,pf] = WaveSeperationWesterhofFunc
    % Perform wave seperation according to the Westerhof et al. 1972
    % method; obtaining pb and pf through calculation of aortic impedance

    global P

    % Initialisation
    Tube          = 'AoSubclAr';
    rho           = 1050; %kg /m3 
    Dt            = P.General.Dt;
    t             = P.t - P.General.tStart;
    tStart        = 1050;
    tOk           = tStart : (tStart + (P.General.tCycle / Dt)); 
    mmHgToPa      = 133.33;

    % calculation [pTrans, c0, Z0] with anti-collaps stiffness of Aortic -
    % Subclavian segment
    Len   = P.Tube.Len   ;                            % repesentative length of blood vessels
    p0    = P.Tube.p0    ;                            % working pressure
    A0    = P.Tube.A0    ;
    rhob  = P.General.rhob;
    A     = max(1e-10,bsxfun(@rdivide,P.Tube.V,Len)); % vessel cross-section
    ANorm = bsxfun(@rdivide,A,A0);                    % cross-section normalized to physiologic volume
    ap    = 0.02;                                     % introduces: Small volume -> negative transmural pressure pTrans
    mk    = 1 + (P.Tube.k/3-2)/(1+ap);                % stiffness exponential
    pp    = (1+ap) * bsxfun(@power, ANorm, mk);
    pm    = -ap./ANorm;                               % anti-collapse pressure
    pTrans= bsxfun(@times,pp + pm, p0);
    c0    = sqrt( bsxfun(@times, bsxfun(@times,pp,mk)-pm, p0) / rhob );
    c0            = c0(:,1); c0 = c0(tOk); %MK wave speed in the aorta
    c0            = rms(c0);
    %

    
    
    % Parsing blood pressure
    AoBP          = Get('Node','p','Ao'); 
    AoBP          = AoBP(tOk);
    AoArea        = Get('Tube','A',Tube); AoArea = AoArea(tOk);
    AoValveQ      = Get('Valve','q','LvAo'); AoValveQ  = AoValveQ(tOk);
    

    % AoQ is sampled differently -> correct for this observation
    AoQ           = Get('Tube','q',Tube); 
    AoQ           = AoQ(tOk);
    AoU           = AoQ / AoArea;
        
    % Calculate characteristic impedance, Zc = p/q [Pa.s /m^3]
    Zc = (rho * c0) ./ AoArea;
    pb = 0.5 * (AoBP - (Zc .* AoQ));
    pf = 0.5 * (AoBP + (Zc .* AoQ));    
    %
    %%
    
 figure;
 plot(AoBP/133,'k','LineWidth',2); hold on;
 plot(pb/133,'b','LineWidth',2); 
 plot(pf/133,'r','LineWidth',2)
 legend('P_s','P_f','P_b')

 P0 = min(AoBP)/133;

figure;
plot(AoBP/133,'c--','LineWidth',2); hold on;
plot( ((pb-pb(1))/133)+((pf-pf(1))/133) + P0 ,'k','LineWidth',1); 
plot(((pb-pb(1))/133),'b','LineWidth',2); 
plot(((pf-pf(1))/133),'r','LineWidth',2)
legend('P_s','P_f + p_b','P_b','P_f')
    
    
end

