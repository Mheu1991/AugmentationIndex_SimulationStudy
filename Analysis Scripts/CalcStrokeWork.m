function LVSW = CalcStrokeWork

% Calculate global stroke work (SW) of the Left Ventricle
%Maarten Heusinkveld 29-06-2015

global P

addpath(genpath('D:\Documents\Universiteit\PhD Jaar 1\Artikel Augmentation Index\CircAdapt\'));


%
iS     = 1050;
dt     = P.General.Dt;
tOk    = iS:(iS + P.General.tCycle/dt);
%
Vlv    = Get('Cavity','V','Lv');
Plv    = Get('Cavity','p','Lv');
%
Vlv    = Vlv(tOk);
Plv    = Plv(tOk);
%

LVSW     = trapz(Plv,Vlv);  %Work in [J/cycle]


end



