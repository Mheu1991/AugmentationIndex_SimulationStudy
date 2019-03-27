%% Obtain parameters from central BP curve and LV ejection and work
% characteristics. Display association between stifness and contractility
% using landscape plots
% Maarten Heusinkveld 16-03-2019

clear all;
close all;
clc;

pathname            = '..\Simulations';

addpath('..\CircAdapt Model')

addpath(pathname);
fileNames           = ls(pathname);
fileNames(1:2,:)    = [];
numFilesToProcess   = size(fileNames,1);
kInit               = 8;
tStart              = 1050;

for i = 1:numFilesToProcess
    load(strtrim(fileNames(i,:)));
    vMax               = Get('Patch','vMax','Lv1');
    kCur               = Get('Tube','k','AoSubclAr');
    dk                 = kCur - kInit;
    
    [PP_ao,PP_ba,Tej,cfPWV,AIx_ao,AIx_Rev_ao,TTP_ao,AIx_ba,AIx_Rev_ba,TTP_ba,RelDn_ao,RelDn_ba,EDV,SV,EDP,LV_SW,AugmentationPressure] = ExtractWaveformProperties(tStart);
    
    CurveProp.(['Sim',num2str(i)]).vMax    = vMax;
    CurveProp.(['Sim',num2str(i)]).dk      = dk;
    CurveProp.(['Sim',num2str(i)]).kCur    = kCur;
    CurveProp.(['Sim',num2str(i)]).PP      = PP_ao;
    CurveProp.(['Sim',num2str(i)]).Tej     = 1e3*Tej;   
    CurveProp.(['Sim',num2str(i)]).cfPWV   = cfPWV;
    CurveProp.(['Sim',num2str(i)]).AIx     = 1e2*AIx_ao; 
    CurveProp.(['Sim',num2str(i)]).AIx_Rev = 1e2*AIx_Rev_ao;
    CurveProp.(['Sim',num2str(i)]).TTP     = 1e3*TTP_ao;  
    CurveProp.(['Sim',num2str(i)]).RelDn   = 1e2*RelDn_ao;
    CurveProp.(['Sim',num2str(i)]).EDV     = EDV;
    CurveProp.(['Sim',num2str(i)]).SV      = SV;
    CurveProp.(['Sim',num2str(i)]).LV_SW   = LV_SW;
    CurveProp.(['Sim',num2str(i)]).AugmentationPressure = AugmentationPressure;

end
   
%%
nofSims = size(fieldnames(CurveProp),1);

v     = zeros(nofSims,1);
k     = zeros(nofSims,1);
aix   = zeros(nofSims,1);
aixrev= zeros(nofSims,1);
tej   = zeros(nofSims,1);
ttp   = zeros(nofSims,1);
lvsw  = zeros(nofSims,1);
augmentationpressure = zeros(nofSims,1);

for i = 1:nofSims   
    v(i)     = CurveProp.(['Sim',num2str(i)]).vMax;
    k(i)     = CurveProp.(['Sim',num2str(i)]).dk;
    aix(i)   = CurveProp.(['Sim',num2str(i)]).AIx;
    aixrev(i)= CurveProp.(['Sim',num2str(i)]).AIx_Rev;
    tej(i)   = CurveProp.(['Sim',num2str(i)]).Tej;
    ttp(i)   = CurveProp.(['Sim',num2str(i)]).TTP;    
    lvsw(i)  = CurveProp.(['Sim',num2str(i)]).LV_SW;
    augmentationpressure(i) = CurveProp.(['Sim',num2str(i)]).AugmentationPressure;
end

dum   = sortrows([v k aix aixrev tej ttp lvsw augmentationpressure],[1,2]);
v     = dum(:,1);
k     = dum(:,2);
aix   = dum(:,3);
aixrev= dum(:,4);
tej   = dum(:,5);
ttp   = dum(:,6);
lvsw  = dum(:,7);
augmentationpressure = dum(:,8);

sK    = size(unique(k),1);
sV    = size(unique(v),1);

AIx     = reshape(aix,sK,sV);
AIx_Rev = reshape(aixrev,sK,sV);
Tej     = reshape(tej,sK,sV);
tFtoS   = reshape(ttp,sK,sV);
Wstroke = reshape(lvsw,sK,sV);
AugmentationPressure = reshape(augmentationpressure,sK,sV);

close all;

figure;
[C,h] = contourf(unique(v),unique(k),AIx);
hold on;
plot(7,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot(7,12,'ks','MarkerSize',10,'MarkerFaceColor','k')
plot(3,0,'ro','MarkerSize',10,'MarkerFaceColor','r')
plot(3,12,'rs','MarkerSize',10,'MarkerFaceColor','r')
xlabel('shortening velocity, v_s [\mu m/s]');
ylabel('\Delta stiffness coefficient, k [-]')
title('AIx [%]')

figure;
[C,h] = contourf(unique(v),unique(k),AIx_Rev);
hold on;
plot(7,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot(7,12,'ks','MarkerSize',10,'MarkerFaceColor','k')
plot(3,0,'ro','MarkerSize',10,'MarkerFaceColor','r')
plot(3,12,'rs','MarkerSize',10,'MarkerFaceColor','r')
xlabel('shortening velocity, v_s [\mu m/s]');
ylabel('\Delta stiffness coefficient, k [-]')
title('AIx Rev [%]')

figure;
[C,h] = contourf(unique(v),unique(k),AIx);
hold on;
plot(7,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot(7,12,'ks','MarkerSize',10,'MarkerFaceColor','k')
plot(3,0,'ro','MarkerSize',10,'MarkerFaceColor','r')
plot(3,12,'rs','MarkerSize',10,'MarkerFaceColor','r')
xlabel('shortening velocity, v_s [\mu m/s]');
ylabel('\Delta stiffness coefficient, k [-]')
title('AIx [%]')

figure;
[C,h] = contourf(unique(v),unique(k),AIx./Tej);
hold on;
plot(7,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot(7,12,'ks','MarkerSize',10,'MarkerFaceColor','k')
plot(3,0,'ro','MarkerSize',10,'MarkerFaceColor','r')
plot(3,12,'rs','MarkerSize',10,'MarkerFaceColor','r')
xlabel('shortening velocity, v_s [\mu m/s]');
ylabel('\Delta stiffness coefficient, k [-]')
title('AIx / t_{ej} [% ms^-1]')

figure;
[C,h] = contourf(unique(v),unique(k),tFtoS);
hold on;
plot(7,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot(7,12,'ks','MarkerSize',10,'MarkerFaceColor','k')
plot(3,0,'ro','MarkerSize',10,'MarkerFaceColor','r')
plot(3,12,'rs','MarkerSize',10,'MarkerFaceColor','r')
xlabel('shortening velocity, v_s [\mu m/s]');
ylabel('\Delta stiffness coefficient, k [-]')
title('time foot-to-shoulder [ms]')

figure;
[C,h] = contourf(unique(v),unique(k),Wstroke);
hold on;
plot(7,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot(7,12,'ks','MarkerSize',10,'MarkerFaceColor','k')
plot(3,0,'ro','MarkerSize',10,'MarkerFaceColor','r')
plot(3,12,'rs','MarkerSize',10,'MarkerFaceColor','r')
xlabel('shortening velocity, v_s [\mu m/s]');
ylabel('\Delta stiffness coefficient, k [-]')
title('LV stroke work [J]')

