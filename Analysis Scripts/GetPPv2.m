function PP = GetPPv2(Pressure)
    % Get PP of aortic BP waveform
    global P
    
    %AObp        = Get('Node','p','Ao')/133;
    tOK         = P.General.tCycle/P.General.Dt;
    PressureTrimmed = Pressure(1:tOK);

    MaxPressure = max(PressureTrimmed);
    MinPressure = min(PressureTrimmed);
    PP          = MaxPressure - MinPressure;
end