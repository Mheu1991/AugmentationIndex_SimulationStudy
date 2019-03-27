function EDP = GetEDP 
    % Estimate EDP of LV 
    % Input is a P-struct of a single beat/ cardiac
    % cycle
    % Maarten Heusinkveld 11-02-2015
    
    global P
    
    dt               = P.General.Dt;
    tOK              = P.General.tCycle/dt;
    LVP              = Get('Cavity','p','Lv');
    
    % obtain a single beat, 2nd order derivative of aortic transvalvular flow
    QAodot                         = Get('Valve','qDot','LvAo')*1e6;
    QAodot                         = QAodot(1:tOK);
    doubleDot                      = diff(QAodot);
    QAodoubleDot(2:length(QAodot)) = doubleDot; 
    QAodoubleDot(1)                = QAodoubleDot(2);

    TstartEjection                 = find( QAodoubleDot == max(QAodoubleDot));
    EDP                            = LVP(TstartEjection);
    
end