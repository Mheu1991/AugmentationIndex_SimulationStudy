function CircAdapt
% function CircAdapt
% Core of the CircAdapt model
% Simulates a series of beats based on de contents of structure P
% Theo Arts, Maastricht University, Oct 13, 2012

global P 
P.General.tStart= P.t(end); % time reference for array counter
P.General.tEnd  = P.t(end)+P.General.DtSimulation;
LinkMatrix; % matrices defining connections between elements

% Initialization P.SVar (= state variables) to last state
P.SVar=[];
P2SVar;
SVar  = P.SVar(end,:);
P.SVar= SVar; % start condition
P.t   = P.t(end);
Dt    = P.General.Dt;
% clearing In-Out record, used for Fast adaptation
P.General.In     = [];
P.General.Out    = []; % reset of storage input/output per beat
P.General.FacpControl= 1; % Initialization pCurrent/pTarget
% if FacpControl>1 => leak of blood volume in ArtVen's

%TimingInit; % Setting depolarization sequence in P.Net.Depolarization
Indexation; % Indexation to define connections between elements

%==== simulation of successive cardiac cycles
while P.t(end)<P.General.tEnd;

    % setting time points for the beat
    tCycle    = P.General.tCycle;
    nt        = round(tCycle/Dt); % number of time points
    TimingInit; % Setting depolarization sequence in P.Net.Depolarization
    
    % Allocation additional storage space for delayed tube signals
    ntExtra=nt+size(P.t,1)-size(P.Tube.pP,1)+2;
    Mat=zeros(ntExtra,P.Tube.n);
    P.Tube.pP = [P.Tube.pP ; Mat];
    P.Tube.pD = [P.Tube.pD ; Mat];
    P.Tube.q  = [P.Tube.q  ; Mat];
    % existing solution
    disp(['t= ',num2str(P.t(end)),';  Time to go= ',...
        num2str(P.General.tEnd-P.t(end))]); pause(0.01);
    % Differential equation
    SVarAppend= OdeCA('SVarDot',Dt,tCycle,P.SVar(end,:));
    SVar= [SVar(1:end-1,:);SVarAppend(1:end,:)]; %appends 1-beat SVar-vector
    % SVar(end,:) removed because of overlap with next beat
    P.SVar= SVarAppend;
    CircDisplay; % display and time course of 1-beat hemodynamics
    % Calculation of changes due to adaptation
    Adapt; % Execute Adapt with P.Adapt.FunctionName
    % Time courses in Par belong to parameter setting before Adapt-action!
end

P.SVar= SVar; %state variables of all beats to show trends. Be careful,
% because AdaptXX has been applied, the P-parameters do not belong to the
% Par.SVar. Only for Adapt0P (= no adaptation), complete SVar
% is compatible with parameter settings.
end

