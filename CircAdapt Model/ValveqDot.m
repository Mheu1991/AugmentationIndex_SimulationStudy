function ValveqDot
% function ValveqDot
% Node pressure differences -> valve flow acceleration qDot
% Valves AvValves refer to ventricular valves with papillary muscles
% Theo Arts, Maastricht University, Jan 1, 2014

global P

iNodeProx=P.Valve.iNodeProx;
iNodeDist=P.Valve.iNodeDist;
nt=size(P.t,1);
AOpen = P.Valve.AOpen; % open valve cross-section
ALeak = P.Valve.ALeak; % closed valve cross-section
AMax  = max(AOpen,ALeak);
%=================
% allows AV-valve diastolic regurgitation
Ws=P.Valve.AvWalls; Vs=P.Valve.AvValves;
T    = P.Wall.T(:,Ws); %wall tension realted to pap. muscle
DADT = P.Wall.DADT(:,Ws); %wall stiffness
Am0  = P.Wall.Am0(:,Ws); %wall zero-stress area
Diast= tanh( 100*max(0.001,T.*DADT./Am0-0.05).^2 ); % diastole->1.0
%=================

% Valve,slow closure to avoid numerical problems
q     = P.Valve.q  ; % flow
Len   = P.Valve.Len; % effective length of flow channel
rhob  = 1050       ; % density of blood
Dp    = P.Node.p(:,iNodeProx)-P.Node.p(:,iNodeDist); % pressure drop
AProx = bsxfun(@minus,P.Node.A(:,iNodeProx),AMax);
ADist = bsxfun(@minus,P.Node.A(:,iNodeDist),AMax);
AMin  = min(AProx,ADist); % minimum external orifice area, fu(t)

AOpenEff = bsxfun(@min,AMin,AOpen); % open flow orifice area, fu(t)
ALeakEff = bsxfun(@min,AMin,ALeak); % leak flow orifice area, fu(t)
ALeakEff(:,Vs)= 0.3*AOpenEff(:,Vs).*Diast+ALeakEff(:,Vs); 
%    Large AV-valve leak in diastole

% Slow closure to avoid numerical problems
Fwq = +(q>0); % Forward flow
Bwq = 1-Fwq ; % Backward flow 
pq  = 50*rhob*q.*abs(q) ./ ( Fwq .* AOpenEff.^2 + Bwq .* ALeakEff.^2 );
%    Flow related pressure term, large factor 50~1/size curved tail
z =(2/pi)*atan2(Dp,pq); % dimensionless state definition of valve
s1= sign(z); % (+/-) represents (forward/backward) flow
y = 1 - 2*sqrt(max(mod(z,2)-1,0)); % square-root steepness at q~0
Fw= s1.*y; % -1<Fw<+1: interpolation factor ALeak <--> AOpen
A = 0.5*((AOpenEff+ALeakEff)+Fw.*(AOpenEff-ALeakEff)); % Valve area

DpB= 0.5*rhob*q.*abs(q).*( A.^-2 - 1./(Fwq.*AProx.^2+Bwq.*ADist.^2) );
%   Bernouilli pressure drop over valve orifice
L  = 1.5*rhob*( bsxfun(@ldivide,A,Len) ...
    + 0.5*(1./sqrt(AProx)+1./sqrt(ADist))); % inertia

P.Valve.L   = L; % inertia
P.Valve.qDot= (Dp-DpB)./L; % flow derivative
P.Valve.AEst = A;

end

