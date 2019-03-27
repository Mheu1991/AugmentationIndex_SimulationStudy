function PP = GetPP
    % Get PP of aortic BP waveform
    global P
    
    AObp        = Get('Node','p','Ao')/133;
    tOK         = P.General.tCycle/P.General.Dt;
    AObpTrimmed = AObp(1:tOK);

    MaxPressure = max(AObpTrimmed);
    MinPressure = min(AObpTrimmed);
    PP          = MaxPressure - MinPressure;
end