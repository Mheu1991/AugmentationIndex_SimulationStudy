function D = WKfuncReservoir(x,x0,t,dt,P0,iNotch,Phat,q)
% Solve equation: q = C dp/dt + p/Rp
% OR            : dp/dt = (q- (p/Rp)/C 
% OR            : iwp   = (q-p/Rp)C
% return pressure in SI units

% C   = 3.5;   %[mL/kPa]   -> x(1)
% Rp  = 0.10;  %[kPa.s/mL] -> x(2)
% Tau = Rp*C;  %[s]

% source data to be fitted
pplot=Phat;
Phat = Phat(iNotch : end); 
%

% transform paraemters back
x = x .* x0;
%

k    = (1+(dt/ (x(2)*x(1))));
%p(1) = P0;
p(1) = min(Phat);

for n = 2:length(t)
    
   p(n) = p(n-1)/k + ( (dt / x(1) ) *q(n) ) /k;
    
end

% Diastolic pressure
ppruned  = p(iNotch:end);

% Vector to be minimized
D = (ppruned' - Phat) / Phat(1);

% plot(t-t(1),p*7.5)
% hold on
% plot(t-t(1),pplot*7.5,'k','LineWidth',2)

end

