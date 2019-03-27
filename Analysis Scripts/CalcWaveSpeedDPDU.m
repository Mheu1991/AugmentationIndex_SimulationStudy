function [c0,P_sys,U_sys,ftk] = CalcWaveSpeedDPDU(Pao,U,MaxI,rho)
% Approximate wave speed from slope of the pressure velocity relation

    P_sys   = Pao(1:MaxI-50);
    U_sys   = U(1:MaxI-50);
    
    s = fitoptions('Method','NonlinearLeastSquares',...
               'Startpoint',[0 2]);
    f = fittype(' a1 + a2 * x ','options',s);
    ftk = fit(U_sys,P_sys,f); % fitting to model P_sys = a1 - a2*U_sys^a3
    
    c0 = ftk.a2 / rho;
    
    %c0 = (P_sys(end) - P_dia(1)) ./ (U_dia(end) - U_dia(1))/rho; 

end

