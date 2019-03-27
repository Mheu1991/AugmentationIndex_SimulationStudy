function Adapt
%function Adapt
% Common to all adaptation procedures.
% Embeds more dedicated adaptation procedures
%    Executed between beats:
% Blood pressure is controlled by change of circulating volume
% Flow is set by P.General.SaturationControl:
%     0-> wanted flow per ArtVen-element
% or: 1-> flow is determined by oxygen usage= AV-Dsaturation x flow
% Theo Arts, Maastricht University, Feb 2, 2014

global P

save PTemp P; %saves last intermediate solution

% Systemic blood flow adjusted per ArtVen-element by resistance change
% Control of systemic pressure by adjustment of circulatory blood volume

%%%%++++%%%

FbFac= P.General.AdaptFeedback; % Feedback constant. If FbFac==0, no control

%%%%++++%%%

% Finding systemic and pulmonary ArtVen's
aux= regexpi(lower(P.ArtVen.Name),'pu'); %all names with 'pu' are pulmonary
Sy = P.ArtVen.Name(cellfun('isempty',aux));% {'Ca','Br','Ce','Fe'} systemic ArtVens
Pu = P.ArtVen.Name(~cellfun('isempty',aux)); % {'Rpu','Lpu'}

%=== Estimation of TargetFlow in systemic ArtVen's for flow control
qSy   = Get('ArtVen','q',Sy);% flows in systemic ArtVen
qPu   = Get('ArtVen','q',Pu);% flows in pulmonary ArtVen
if P.General.SaturationControl
    % 0: flow set / 1: control of flow by AV-saturation difference
    % Preparation execution of Saturation.m
    % Setting venous target saturation allowing execution of Saturation.m
    % Quantitative specific data in program text ++++++++++
    SatVenSy=[0.50,0.75,0.75,0.75,0.75,0.75,0.75,0.75]; %{'Ca','Br','Ce','Fe'} +++++++++{'Ca','Br1','Br2','Br3','Fi','Ce','Fe'};
    SatVenPu=0.98;
    Put('ArtVen','SatVen',Sy,SatVenSy);
    Put('ArtVen','SatVen',Pu,SatVenPu);
    P.General.nBeatSat=5; % number of cycles for to approach steady state of saturation
    Saturation
    
    % TargetFlow qSyTarget determined by systemic oxygen consumption
    % expressed as qDSatSyTarget= flow*DSaturation(A-V) in systemic
    % ArtVen's
    qDSatSyTarget  = [0.072,0.032,0.131,0.037]*P.General.q0;
    SatSyProx      = P.Node.Sat(:,Get('ArtVen','iNode',Sy)); % proximal node saturation
    qDSatSy        = mean((SatSyProx-ones(length(SatSyProx),1)*SatVenSy).*qSy); % mean flow*DavSaturation product
    FacqDSatControl= exp((1-qDSatSy./qDSatSyTarget)); % Vasodilation factor
    qSyTarget      = Get('ArtVen','q0AV',Sy).*FacqDSatControl; % Target flow in systemic ArtVen's
    % Vasodilation if consumption > delivery by actual flow
else % direct flow control
    %Armq0 = 0.12;
    qSyTarget=[0.15,0.01,0.07,0.02,0.015,0.005,0.57,0.16]*P.General.q0; %+++++++ Target flow in systemic ArtVen's [15,12,57,16] [0.15,Armq0/3,Armq0/4,Armq0/6,15*Armq0/64,Armq0/64,0.49,0.16]
    %qSyTarget = 1e-2*[15,12,57,16]*P.General.q0; 
end

q0AV=Get('ArtVen','q0AV',Sy);
Put('ArtVen','q0AV',Sy,q0AV.*(qSyTarget./q0AV).^FbFac); % save Target flows in P.ArtVen.q0AV
%=== End of Estimation of TargetFlow

% Control of flow to P.ArtVen.q0AV for systemic ArtVen's
q0AvSy=Get('ArtVen','q0AV',Sy);
p0AvSy=Get('ArtVen','p0AV',Sy);
kAvSy =Get('ArtVen','kAV' ,Sy);
FacqControlSy = exp(FbFac*(1-mean(qSy)./q0AvSy));
Put('ArtVen','p0AV',Sy,p0AvSy.*FacqControlSy.^(-1./kAvSy));

% Flow in Sys and Pu to judge steady state
FlowVec=[sum(mean(qSy)),sum(mean(qPu))]; % Systemic/Pu flow
disp('Flow/q0: Sys Pu');
disp(FlowVec/P.General.q0);

%Control of pressure by change of circulating volume
FacpControl=mean(P.Node.p(:,P.Node.iBaro))/P.General.p0;
P.General.FacpControl=FacpControl^FbFac;

% Estimate AV-delay
P.General.TauAv=0.185*P.General.tCycle;

% ====== AdaptSpecific, different for Rest, Exercise
NoAdapt= strcmp(P.General.AdaptFunction,'Adapt0P');
if NoAdapt % No adaptation, regular sequence of beats 
    % === Faster Steady State at rest
    VecV=Get('Cavity','V','All');
    P.General.In =[P.General.In ;VecV(  1,:) ];
    P.General.Out=[P.General.Out;VecV(end,:) ];
    if P.General.Fast;% Put new start values in structure P
        Vec= SteadyStateP;
        i=1  ; j=P.Cavity.n; Put('Cavity','V','All',Vec(i:j));
    end
else
    feval(P.General.AdaptFunction); % specific adapt function
end

% Judging quality of steady state
if size(P.General.Out,1)>1; % Escape if steady state is reached
    ErrVec= 1000*log( P.General.Out(end,:)./P.General.In(end,:) );
    disp(['Stationarity error= ',...
        num2str(round(norm(ErrVec .^2 )))] );
    %=== ERROR criterium on flow stationarity
    if sqrt(mean(ErrVec .^ 2 ))<1 && (~NoAdapt || P.General.Fast)
        P.General.tEnd= P.General.tEnd-0.5*(P.General.tEnd-P.t(end));
    end
end;

% get the initial condition for next beat, P is most compact information to
% start the simulation
P2SVar; % load physiologic data in record of state variables P.SVar

end

