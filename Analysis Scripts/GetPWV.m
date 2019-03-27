function cfPWV = GetPWV

    %% Acquire PWV in tubes and/ or ArtVen element per simulation
    % Maarten Heusinkveld 4-06-2014
    %%
    
    global P
    
    Seg    = {'CaAr','FeAr'};
    SegLen = 0.43;
    dt     = P.General.Dt;
    samp   = 1/dt;
    onsetT = zeros(1,size(Seg,2));

    if length(Get('Node','p','Ao')) == 827
        iS  = 130;
    else
        iS   = 1080;
    end
    
    tOK  = iS:(iS + P.General.tCycle/dt);

    
    for i = 1:size(Seg,2)

        %Calc 2nd derivative
        datS = Get('Node','p',Seg{i});    
        d1S = diff(datS)/dt;
        D1S(2:length(datS)) = d1S; D1S(1) = D1S(2);
        d2S = diff(D1S)/dt; 
        D2S(2:length(D1S))  = d2S; D2S(1) = D2S(2);

        onsetT(i) = find(D2S == max(D2S(tOK)));

    end

    cfPWV = (SegLen / ( (onsetT(2) - onsetT(1)) /samp ));

    %figure;
    Ca     = Get('Node','p',Seg{1})/133;
    Fe     = Get('Node','p',Seg{2})/133;

%     plot(Ca); hold on; plot(Fe,'r')
%     plot(onsetT(1),Ca(onsetT(1)),'ks'); plot(onsetT(2),Fe(onsetT(2)),'rs')
% 
%     disp(['cfPWV = ',num2str(cfPWV) ' m/s'])
end
