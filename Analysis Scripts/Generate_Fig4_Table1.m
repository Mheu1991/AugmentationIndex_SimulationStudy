%% Make table summarizing simulation outcome for the follwoing simulation scenarios
% - Reference model
% - Vessel stiffening
% - Reduced LV contractility
% - Combined

% Ps, Pd, Pp, Pm : systolic, diastolic, pulse, mean pressure
% RelDn_Ao, RelDn_Ba: relative height of the aortic, brachial dicrotic
% notch
% I_aAo, I_aBa, I_Rev_aAo: aortic augmentation index, brachial augmentation
% indec, alternative augmentation index
% Tej: LV ejection duration
% TTP: time-to-inflection
% cfPWV: carotid-femoral pulse wave velocity
clear all;
close all;
clc;

addpath(genpath('..\Simulations\'));

FileN    = {'PNewvMax7kTube0TubeArtVen.mat','PNewvMax7kTube12TubeArtVen.mat','PNewvMax3kTube0TubeArtVen.mat','PNewvMax3kTube12TubeArtVen.mat'}; % <-- manuscript sims 

NumFiles = size(FileN ,2);
mmHgToPa = 133.33;
tStart   = 1050;

[Ps,Pd,Pp,Pm]                          = deal(zeros(NumFiles,1));
[RelDn_Ao,RelDn_Ba]                    = deal(zeros(NumFiles,1));
[I_aAo,I_aBa,I_Rev_aAo]                = deal(zeros(NumFiles,1));
[Tej,TTP]                              = deal(zeros(NumFiles,1));
[EDV,SV]                               = deal(zeros(NumFiles,1));
[pbpf_ws,pbpfpb_ws]                    = deal(zeros(NumFiles,1));
[pb_ws,pf_ws]                          = deal(zeros(NumFiles,1));
[pbpf_fou,pb_fou,pf_fou]               = deal(zeros(NumFiles,1));
LVSW                                   = zeros(NumFiles,1);
PSSR                                   = zeros(NumFiles,1);
   
for i = 1:NumFiles
load(FileN{i})
Ps(i)                                                                                    = max(Get('Node','p','Ao'))/mmHgToPa;
Pd(i)                                                                                    = min(Get('Node','p','Ao'))/mmHgToPa;
Pm(i)                                                                                    = mean(Get('Node','p','Ao'))/mmHgToPa;
[Pp(i),~,Tej(i),~,I_aAo(i),I_Rev_aAo(i),TTP(i),~,~,~,~,~,EDV(i),SV(i),~,LVSW(i)]         = ExtractWaveformProperties(tStart);
[pbpf_ws(i),pbpfpb_ws(i),~,~,~,pb_ws(i),pf_ws(i),~,~]                                    = PressureWaveSeparationFunction('wavesep','true',tStart);
PSSR(i)                                                                                  = StrainRateFunc;

end

Tej      = 1e3* Tej;
I_aAo    = 1e2 * I_aAo;
I_aBa    = 1e2 * I_aBa;
RelDn_Ao = 1e2* RelDn_Ao;
RelDn_Ba = 1e2* RelDn_Ba;
I_Rev_aAo= 1e2* I_Rev_aAo;
EDV      = 1e6* EDV;
TTP      = 1e3 * TTP;

sName = {'Reference','Stiffening','Reduced vMax','Combined'};
vName = {'p_sys','p_dia','p_pulse','p_mean','I_a','I_Rev_aAo','pbpf_WS','pb_WS','pf_WS','TTP','Tej','LVSW','PSSR'};
dat   = [Ps,Pd,Pp,Pm,I_aAo,I_Rev_aAo,pbpf_ws,pb_ws,pf_ws,TTP,Tej,LVSW,PSSR]';
        
table(Ps,Pd,Pp,Pm,I_aAo,I_Rev_aAo,pbpf_ws,pb_ws,pf_ws,TTP,Tej,LVSW,PSSR,'VariableNames',vName,'RowNames',sName)



