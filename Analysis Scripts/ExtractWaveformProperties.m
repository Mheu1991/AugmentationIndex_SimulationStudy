function [PP_ao,PP_ba,Tej,cfPWV,AIx_ao,AIx_Rev_ao,TTP_ao,AIx_ba,AIx_Rev_ba,TTP_ba,RelDn_ao,RelDn_ba,EDV,SV,EDP,LV_SW,AugmentationPressure] = ExtractWaveformProperties(tStart)
    
    global P

    %Initialisation
    mmHgToPa= 133.33;
    %tStart  = 1050;
    tOk     = tStart : (tStart + (P.General.tCycle / P.General.Dt)); 
    AoBP    = Get('Node','p','Ao')/mmHgToPa;
    BraBP    = Get('Node','p','BraAr')/mmHgToPa;
    %AoBP    = AoBP(tOk);
    QAoDot  = Get('Valve','qDot','LvAo')*1e6;
    QAoDot  = QAoDot(tOk);
    %
    PP_ao   = GetPPv2(AoBP);
    PP_ba   = GetPPv2(BraBP);
    Tej     = EjectionTime(QAoDot);
    cfPWV   = GetPWVv2;
    [AIx_ao,TTP_ao,AIx_Rev_ao,~,AugmentationPressure]  = AugmentationIndex_v4(AoBP);
    [AIx_ba,TTP_ba,AIx_Rev_ba]  = AugmentationIndex_v4(BraBP);
    RelDn_ao= RelDnHeight(AoBP);
    RelDn_ba= RelDnHeight(BraBP);
    EDV     = GetEDV;
    SV      = GetSV;
    EDP     = GetEDP;
    LV_SW   = CalcStrokeWork;
end