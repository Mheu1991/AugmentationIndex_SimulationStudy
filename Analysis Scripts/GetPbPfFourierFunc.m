function [pbpf,Rc ] = GetPbPfFourierFunc( Analysis,PlotYoN,tStart )

% Seperate waveform in forward and backward component according to
% Westerhof et al. 1972 )

% Initialisation

global P

dt          = P.General.Dt;
sf          = 1/dt;
tStart      = 1;
mmHgToPa    = 133.33;
freq        = 1/(P.General.tCycle); %[Hz]
noHarmonics = 10;
fComp       = [0:1:noHarmonics] * freq;
ns          = 1e3;

% define time vector of simulated data (tS) 
tS          = P.t - P.General.tStart;

% parse simulated pressure and flow
AoBPs       = Get('Tube','p','AoSubclAr');
pd          = min(AoBPs);
AoBPs       = AoBPs - pd;
AoBP        = AoBPs;
AoQs        = Get('Tube','q','AoSubclAr');  
AoQ         = AoQs;

Pcoef      = fourier(AoBPs,sf,noHarmonics,freq);
Qcoef      = fourier(AoQs,sf,noHarmonics,freq);
Zcoef      = Pcoef ./ Qcoef;
ZcoefAngle = angle(Zcoef);

%%
% Seperate waveform in forward and backward component according to
% Westerhof et al. 1972 )

% Zc: characteristic impedance is taken as the minimal of Zin in the 0-10
% Hz interval

Zc             = mean(Zcoef(end-5:end));

% calculate the moduli of pf and pb
for n = 1: noHarmonics+1
    
    Z          = Zcoef(n);        
    
    Rc(n)      = (Z - Zc) / (Z + Zc);    
          
    PFcoef(n)      = Pcoef(n) /(1+Rc(n));
    PBcoef(n)      = (Pcoef(n) * Rc(n)) / (1+Rc(n));
    
    QFcoef(n)      = Qcoef(n) /(1-Rc(n));
    QBcoef(n)      = -Rc(n) * (Qcoef(n) / (1-Rc(n))); 
          
end


%% inverse fourier transformation to time-domain signal

pf = ifourier(PFcoef, ns, freq);
pb = ifourier(PBcoef, ns, freq);
qf = ifourier(QFcoef, ns, freq);
qb = ifourier(QBcoef, ns, freq);

%% Plotting

figure;
subplot(2,1,1); %hold on;
plot(fComp,Zcoef,'ks--','LineWidth',1.5)
xlabel('frequency [Hz]')
ylabel('|Z|')
axis('tight')

%figure;
subplot(2,1,2);
plot(fComp,ZcoefAngle,'ks--','LineWidth',1.5);
hold on;
xlabel('frequency [Hz]')
ylabel('phase [rad]')
axis('tight')
legend('Reference')

figure;
plot(pf/mmHgToPa,'b'); hold on;
plot(pb/mmHgToPa,'r');
plot((pf+pb)/mmHgToPa,'k')
legend('pf','pb','pf+pb','ps')

figure;
plot(qf,'b'); hold on;
plot(qb,'r');
plot(qf+qb,'k')
legend('qf','qb','qf+qb','qs')

figure;
subplot(4,1,1)
stem(fComp,Rc)
hold on;
xlabel('frequency [Hz]')
ylabel('RC [-]')

subplot(4,1,2)
stem(fComp,PFcoef/mmHgToPa)
xlabel('frequency [Hz]')
ylabel('PF [mmHg]')

subplot(4,1,3)
stem(fComp,PBcoef/mmHgToPa)
xlabel('frequency [Hz]')
ylabel('PB [mmHg]')

subplot(4,1,4)
stem(fComp,(PBcoef+PFcoef)/mmHgToPa)
xlabel('frequency [Hz]')
ylabel('PB + PF [mmHg]')

%% calculate relection index amplitude(pb)/amplitude(pf)

pbA = max(pb)/ mmHgToPa
pfA = max(pf)/ mmHgToPa

pbpf = pb/pf

end

