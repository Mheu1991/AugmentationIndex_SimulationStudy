function Tej = EjectionTime(QAodot) 
    % Estimate time of ejection from LV in the aorta based on flow profiles
    % over the aortic valve. Input is a P-struct of a single beat/ cardiac
    % cycle
    % Maarten Heusinkveld 27-06-2014
    
    global P
    
    dt               = P.General.Dt;
    tOK              = P.General.tCycle/dt;
    RestEjectionTime = 0.300; % normal ejection time ~300 ms at rest ; find second maximum in qDoubleDot for end ejection

    % obtain a single beat, 2nd order derivative of aortic transvalvular flow
    %QAodot                         = Get('Valve','qDot','LvAo')*1e6;
    %QAodot                         = QAodot(1:tOK);
    doubleDot                      = diff(QAodot);
    QAodoubleDot(2:length(QAodot)) = doubleDot; 
    QAodoubleDot(1)                = QAodoubleDot(2);

    IstartEjection                 = find( QAodoubleDot == max(QAodoubleDot));
    MidEjEstimated                 = (0.5*RestEjectionTime) / P.General.Dt; 
    SecPart                        = MidEjEstimated + IstartEjection;

    IstopEjection                  = find( QAodoubleDot(SecPart:end) == max(QAodoubleDot(SecPart:end))) + SecPart;

    Tej                            = (IstopEjection -  IstartEjection) * P.General.Dt;
    
%     figure; plot(QAodoubleDot); hold on;
%     plot(IstartEjection,QAodoubleDot(IstartEjection),'ks'); 
%     plot(IstopEjection,QAodoubleDot(IstopEjection),'ks')
end


