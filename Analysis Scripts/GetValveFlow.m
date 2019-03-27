addpath('D:\Documents\Universiteit\Jaar 5\Afstuderen\CircAdaptO212a\Simulations Artery14');
load('PVmax4_Aging_6PressureVolumeOff.mat');

dT       = P.General.Dt;
Aoq      = Get('Valve','q','LvAo')*1e6;

plot(Aoq);
xlabel('Iteration [-]')
ylabel('Aortic valve flow [mL/s]')

[x,~]    = ginput();
X        = round(x);
AoqBeat  = Aoq(X(1):X(2));

AortaOut = dT * sum(AoqBeat)