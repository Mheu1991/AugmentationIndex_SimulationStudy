function pRes = WKfunc(C,Rp,t,dt,P0,q)
% Solve equation: q = C dp/dt + p/Rp
% OR            : dp/dt = (q- (p/Rp)/C 
% OR            : iwp   = (q-p/Rp)C
% return pressure in SI units

% C   = 5.3483;   %[mL/kPa]
% Rp  = 0.1423;  %[kPa.s/mL]

k    = (1+(dt/ (Rp*C )));

p(1) = P0;

for n = 2:length(t)
    
   p(n) = p(n-1)/k + ( (dt/C) *q(n) ) /k;
    
end

pRes = p *1e3;

end

