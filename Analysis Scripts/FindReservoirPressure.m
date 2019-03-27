function [CFitted,RpFitted,ResP,P0F,resPerc] = FindReservoirPressure(AoBP,AoBP1cycle,AoQ1cycle,t1cycle,Dt)
    
    global P
    
    % Fit reservoir pressure to CircAdapt simulated central BP curve
    % LSQ method

    % find dicrotic notch time and pressure (~ valve closure)
    [~ , ~, DnPos] = RelDnHeight(AoBP1cycle);
    P0F            = min(AoBP) / 1e3;   % kPa
    AoBP1cycleF    = AoBP1cycle / 1e3 ; % kPa
    AoQ1cycleF     = AoQ1cycle * 1e6;   % mL
    %

    % initialisation fitting
    
    
    x0 = [3.8, (1e-3*P.General.p0) / (1e6*P.General.q0)];
    lb = [0.01, 0.01];
    ub = [50, 2.0];

    % normalisation
    x0N = x0 ./ x0;
    lbN = lb ./ x0;
    ubN = ub ./ x0;

    % LSQ fitting
    opts                  = optimoptions(@lsqnonlin,'Display','off');
    [xN,resnorm,residual] = lsqnonlin(@(x)WKfuncReservoir(x,x0,t1cycle,Dt,P0F,DnPos,AoBP1cycleF,AoQ1cycleF),x0N,lbN,ubN,opts);
    

    resPerc = resnorm / P0F;

    % Back to fitted parameters
    x    = xN .* x0;

    CFitted = x(1); %[mL/kPa]
    RpFitted= x(2); %[kPa.s/mL]

    ResP = WKfunc(x(1),x(2),t1cycle,Dt,P0F,AoQ1cycleF);
%
end

