%% Calculate and plot fibre strain rate 
% M Heusinkveld 21-08-2018
clear all; close all; clc;

load('D:\Documents\Universiteit\PhD Jaar 1\20151004 Artikel Augmentation Index\CircAdapt\SimulationsTubeArtVen\PNewvMax7kTube12TubeArtVen.mat')

t      = P.t-P.t(1);
Efiber = exp((P.Patch.Ef(:,4)+P.Patch.Ef(:,4))/2); % logarithmic to cauchy strain
SRfiber= gradient(Efiber')/P.General.Dt;

figure(1);
plot(t,SRfiber,'b-','LineWidth',2); hold on;

Ref_PeakSR = abs(min(SRfiber))

figure(2);
plot(t,Efiber,'b-','LineWidth',2); hold on;

load('D:\Documents\Universiteit\PhD Jaar 1\20151004 Artikel Augmentation Index\CircAdapt\SimulationsTubeArtVen\PNewvMax3kTube12TubeArtVen.mat')

t      = P.t-P.t(1);
Efiber = exp((P.Patch.Ef(:,4)+P.Patch.Ef(:,4))/2); % logarithmic to cauchy strain
SRfiber= gradient(Efiber')/P.General.Dt;

Reduced_PeakSR = abs(min(SRfiber))

figure(1);
plot(t,SRfiber,'r--','LineWidth',2); hold on;
grid on;
set(gca,'FontSize',14)
xlabel('time [s]')
ylabel('strain rate [1/s]')
legend('NORM','Reduced vMax','Location','Best')

figure(2);
plot(t,Efiber,'r--','LineWidth',2); hold on;
grid on;
set(gca,'FontSize',14)
xlabel('time [s]')
ylabel('fiber strain [-]')
legend('NORM','Reduced vMax','Location','Best')

1e2*((Reduced_PeakSR-Ref_PeakSR)/Ref_PeakSR)
