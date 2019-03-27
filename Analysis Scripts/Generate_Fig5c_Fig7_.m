%% Extract pressure wavefrom charactersistics from simulation batch and visualize in contour (landscape-) plot
clear all;
close all;
clc;

%%
pathname            = '..\Simulations';
addpath(pathname);
addpath('..\CircAdapt Model\');
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
    
    [pbpf,pbpfpb,CFitted,RpFitted,ErrFitted,pb,pf,C0,pres_integrated] = PressureWaveSeparationFunction('wavesep','false',tStart);
    [WRI,SwaveArea,CwaveArea]  = WaveReflectionIndexFunc('false',tStart);
    
    CurveProp.(['Sim',num2str(i)]).vMax            = vMax;
    CurveProp.(['Sim',num2str(i)]).dk              = dk;
    CurveProp.(['Sim',num2str(i)]).kCur            = kCur;
    CurveProp.(['Sim',num2str(i)]).pbpf            = pbpf;
    CurveProp.(['Sim',num2str(i)]).pbpfpb          = pbpfpb;
    CurveProp.(['Sim',num2str(i)]).pb              = pb;
    CurveProp.(['Sim',num2str(i)]).pf              = pf;
    CurveProp.(['Sim',num2str(i)]).WRI             = WRI;
    CurveProp.(['Sim',num2str(i)]).CwaveArea       = CwaveArea;
    CurveProp.(['Sim',num2str(i)]).SwaveArea       = SwaveArea;
    CurveProp.(['Sim',num2str(i)]).CFitted         = CFitted;
    CurveProp.(['Sim',num2str(i)]).RpFitted        = RpFitted;
    CurveProp.(['Sim',num2str(i)]).ErrFitted       = ErrFitted;
    CurveProp.(['Sim',num2str(i)]).C0              = C0;
    CurveProp.(['Sim',num2str(i)]).pres_integrated = pres_integrated;

end

%%
nofSims = size(fieldnames(CurveProp),1);

v                = zeros(nofSims,1);
k                = zeros(nofSims,1);
pbpf             = zeros(nofSims,1);
pbpfpb           = zeros(nofSims,1);
pb               = zeros(nofSims,1);
pf               = zeros(nofSims,1);
wri              = zeros(nofSims,1);
cwave            = zeros(nofSims,1);
swave            = zeros(nofSims,1);
cfit             = zeros(nofSims,1);
rpfit            = zeros(nofSims,1);
errfit           = zeros(nofSims,1);
c0               = zeros(nofSims,1);
pres_integrated  = zeros(nofSims,1);

for i = 1:nofSims   
    v(i)     = CurveProp.(['Sim',num2str(i)]).vMax;
    k(i)     = CurveProp.(['Sim',num2str(i)]).dk;
    pbpf(i)  = CurveProp.(['Sim',num2str(i)]).pbpf;
    pbpfpb(i)= CurveProp.(['Sim',num2str(i)]).pbpfpb;
    pb(i)    = CurveProp.(['Sim',num2str(i)]).pb;
    pf(i)    = CurveProp.(['Sim',num2str(i)]).pf;
  
end

dum   = sortrows([v k pbpf pbpfpb pb pf],[1,2]);
v     = dum(:,1);
k     = dum(:,2);
pbpf  = dum(:,3);
pbpfpb= dum(:,4);
pb    = dum(:,5);
pf    = dum(:,6);


sK    = size(unique(k),1);
sV    = size(unique(v),1);

PBPF  = reshape(pbpf,sK,sV);
PBPFPB= reshape(pbpfpb,sK,sV);
PB    = reshape(pb,sK,sV);
PF    = reshape(pf,sK,sV);

close all;

figure(1);
[C,h] = contourf(unique(v),unique(k),1e2*PBPF);
hold on;
plot(7,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot(7,12,'ks','MarkerSize',10,'MarkerFaceColor','k')
plot(3,0,'ro','MarkerSize',10,'MarkerFaceColor','r')
plot(3,12,'rs','MarkerSize',10,'MarkerFaceColor','r')
xlabel('shortening velocity, v_s [\mu m/s]');
ylabel('\Delta stiffness coefficient, k [-]')
title('p_b/p_f [%]')

figure(2);
[C,h] = contourf(unique(v),unique(k),1e2*PBPFPB);
hold on;
plot(7,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot(7,12,'ks','MarkerSize',10,'MarkerFaceColor','k')
plot(3,0,'ro','MarkerSize',10,'MarkerFaceColor','r')
plot(3,12,'rs','MarkerSize',10,'MarkerFaceColor','r')
xlabel('shortening velocity, v_s [\mu m/s]');
ylabel('\Delta stiffness coefficient, k [-]')
title('p_b/(p_f + p_b)[%]')

figure(9);
[C,h] = contourf(unique(v),unique(k),PB);
hold on;
plot(7,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot(7,12,'ks','MarkerSize',10,'MarkerFaceColor','k')
plot(3,0,'ro','MarkerSize',10,'MarkerFaceColor','r')
plot(3,12,'rs','MarkerSize',10,'MarkerFaceColor','r')
xlabel('shortening velocity, v_s [\mu m/s]');
ylabel('\Delta stiffness coefficient, k [-]')
title('p_b [mmHg]')

figure(10);
[C,h] = contourf(unique(v),unique(k),PF);
hold on;
plot(7,0,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot(7,12,'ks','MarkerSize',10,'MarkerFaceColor','k')
plot(3,0,'ro','MarkerSize',10,'MarkerFaceColor','r')
plot(3,12,'rs','MarkerSize',10,'MarkerFaceColor','r')
xlabel('shortening velocity, v_s [\mu m/s]');
ylabel('\Delta stiffness coefficient, k [-]')
title('p_f [mmHg]')

%%