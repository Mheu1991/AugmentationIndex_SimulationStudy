function cfPWV = GetPWVv2

    %% Acquire PWV in tubes and/ or ArtVen element per simulation
    % Maarten Heusinkveld 6-02-2015
    %%    
    global P
    
    Seg    = {'CaAr','FeAr'};
    SegLen = (Get('Tube','Len','AoSubclAr') + Get('Tube','Len','SubclArCeAr') + Get('Tube','Len','CeArFeAr')) - (Get('Tube','Len','SubclArCaAr') + Get('Tube','Len','AoSubclAr') );
    tCycle = P.General.tCycle;
    onsetT = zeros(1,size(Seg,2));
   
    tStart  = 1050;
    dt      = 1/P.General.Dt;
    tOk     = tStart : (tStart + (P.General.tCycle / P.General.Dt)); 
    
    for i = 1:size(Seg,2)

        %Calc 2nd derivative
        datS = Get('Node','p',Seg{i});    
        datS = datS(tOk);
        d1S = diff(datS)/dt;
        D1S(2:length(datS)) = d1S; D1S(1) = D1S(2);
        eval(['n1S',num2str(i),'=','D1S;']);
        d2S = diff(D1S)/dt; 
        D2S(2:length(D1S))  = d2S; D2S(1) = D2S(2);
        eval(['n2S',num2str(i),'=','D2S;']);
        onsetT(i) = find(D2S == max(D2S(1:(1e3*tCycle))));

    end

    cfPWV = (SegLen / ( (onsetT(2) - onsetT(1)) /dt ));

    %subplot(3,1,1) %figure;
    Ca     = Get('Node','p',Seg{1})/133;
    Ca     = Ca(tOk);
    Fe     = Get('Node','p',Seg{2})/133;
    Fe     = Fe(tOk);

%     subplot(3,1,1)
%     plot(Ca); hold on; plot(Fe,'r')
%     plot(onsetT(1),Ca(onsetT(1)),'ks');
%     plot(onsetT(2),Fe(onsetT(2)),'rs')
%     
%     subplot(3,1,2)
%     plot(n1S1);hold on; plot(n1S2);
%     plot(onsetT(1),n1S1(onsetT(1)),'ks')
%     plot(onsetT(2),n1S2(onsetT(2)),'rs')
%     
%     subplot(3,1,3)
%     plot(n2S1);hold on; plot(n2S2);
%     plot(onsetT(1),n2S1(onsetT(1)),'ks')
%     plot(onsetT(2),n2S2(onsetT(2)),'rs')
%     
%     disp(['cfPWV = ',num2str(cfPWV) ' m/s'])
end
