function Indexation
% function Indexation
% Sets backbone of structure P
% Attribute names to ArtVen, TriSeg, Chambers, Valves, Tubes, 
% Walls and Patches
% Connections between elements are defined by strings stored in
% structure field P.Net
%
% Structure is built on information in P.Net
% ArtVen represents artery-peripheral resistance-vein of organ or body part
% Chamber represents a cavity enclosed by a myocardial wall (atria)
% TriSeg represents combination of two cavities with three walls
% Bag represents passive elastic bag, encapsulating part of the
% circulation, like the pericardium
% Node: named connection point
% Cavity: volume with elastic or muscular wall, connected to a node
% Wall: muscular wall can contract, composed of 1 or more patches
% Patch: contracting part of a wall, having specific mechanical properties
% Valve: valve with inertia connects proximal to distal node, may leak
% Tube: elastic tube connects proximal to distal node
%
% Theo Arts, Maastricht University, April 8, 2013

global P;

Net=P.Net; % contains complete information about element composition 
% and mutual connections of Heart and Circulation

%==== All about Name-strings of elements =============

% ArtVen->Cavity
ArtVenName   = Net.ArtVen.Name;
ArtVen2Cavity= strcat(ArtVenName,{'Ar'});
Aux=[strcat(ArtVenName,{'Ar'});strcat(ArtVenName,{'Ve'})];
CavityName= Aux(:)';

% Chamber->[Cavity,Wall]
ChamberName    = Net.Chamber.Name;
Chamber2Cavity = ChamberName;
CavityName     = [CavityName,Chamber2Cavity];
Chamber2Wall   = Chamber2Cavity;
WallName       = Chamber2Wall;

% TriSeg->[Cavity,Wall]
TriSegName   = Net.TriSeg.Name;
TriSeg2Cavity= strcat({'L'},TriSegName);
Aux= [strcat({'L'},TriSegName);strcat({'R'},TriSegName)];
CavityName   = [CavityName,Aux(:)'];
TriSeg2Wall  = strcat({'L'},TriSegName);
Aux=[strcat({'L'},TriSegName);...
    strcat({'S'},TriSegName);...
    strcat({'R'},TriSegName)];
WallName     = [WallName,Aux(:)'];

% Cavity->Node
Cavity2Node  = CavityName;
NodeName     = Cavity2Node;

% Valve->NodeProx/Dist
ValveNode     = Net.Valve.Nodes;
Valve2NodeProx= ValveNode(:,1)';
Valve2NodeDist= ValveNode(:,2)';
ValveName     = strcat(Valve2NodeProx,Valve2NodeDist);

% Tube->NodeProx/Dist
TubeNode     = Net.Tube.Nodes;
Tube2NodeProx= TubeNode(:,1)';
Tube2NodeDist= TubeNode(:,2)';
TubeName = strcat(Tube2NodeProx,Tube2NodeDist);
dNodeName= setdiff([Tube2NodeProx,Tube2NodeDist],NodeName); % add non-existing node(s)
NodeName = [NodeName,dNodeName];

% Wall->Patch, Wall.nPatch
Wall2Patch= strcat(WallName,'1');
PatchName = Wall2Patch;
WallnPatch=ones(1,length(WallName));
% MultiPatch
MultiPatchWall= Net.Wall.MultiPatch;
nMultiPatch   = length(MultiPatchWall);
iW=(1:length(WallName))';
iP=(1:length(PatchName))';
for i=1:nMultiPatch
    NameMulti={};
    Str=MultiPatchWall{i};
    LtrStr   =isletter(Str);%logical letter positions
    NameStr  =Str(LtrStr);
    NumberStr=str2double(Str(~LtrStr));
    iWall    =strcmp( NameStr     ,WallName )*iW;
    iPatch   =strcmp([NameStr,'1'],PatchName)*iP;
    for k=2:NumberStr
        NameMulti=[NameMulti,[NameStr,num2str(k)]];
    end
    PatchName=[PatchName(1:iPatch),NameMulti,PatchName(iPatch+1:end)];
    WallnPatch(iWall)=NumberStr;
end

BagName    = Net.Bag.Name;

%==== END All about Name strings =============

P.Node.Name=[]; % clears P.Node in DataWrite

%==== indexations ==========
% Write data in columns, determined by name strings
P.ArtVen = DataWrite('ArtVen' ,ArtVenName );
P.Chamber= DataWrite('Chamber',ChamberName);
P.TriSeg = DataWrite('TriSeg' ,TriSegName );
P.Valve  = DataWrite('Valve'  ,ValveName  );
P.Tube   = DataWrite('Tube'   ,TubeName   );
P.Node   = DataWrite('Node'   ,NodeName   );
P.Cavity = DataWrite('Cavity' ,CavityName );
P.Wall   = DataWrite('Wall'   ,WallName   );
P.Patch  = DataWrite('Patch'  ,PatchName  );

% indices determine mutual relations between elements
P.ArtVen.iCavity = Get('Cavity','Index',ArtVen2Cavity );
P.Chamber.iCavity= Get('Cavity','Index',Chamber2Cavity);
P.Chamber.iWall  = Get('Wall'  ,'Index',Chamber2Wall  );
P.TriSeg.iCavity = Get('Cavity','Index',TriSeg2Cavity );
P.TriSeg.iWall   = Get('Wall'  ,'Index',TriSeg2Wall   );
P.Valve.iNodeProx= Get('Node'  ,'Index',Valve2NodeProx);
P.Valve.iNodeDist= Get('Node'  ,'Index',Valve2NodeDist);
P.Tube.iNodeProx = Get('Node'  ,'Index',Tube2NodeProx );
P.Tube.iNodeDist = Get('Node'  ,'Index',Tube2NodeDist );
P.Cavity.iNode   = Get('Node'  ,'Index',Cavity2Node   );
P.Wall.nPatch    = WallnPatch;
P.Wall.iPatch    = Get('Patch' ,'Index',Wall2Patch    );

% additional through-indexation (for convenience)
P.ArtVen.iNode  = P.Cavity.iNode(P.ArtVen.iCavity  );

% identify element content of bags by indices
P.Bag.Name   = BagName;
nBag         = length(P.Bag.Name);
P.Bag.n      = nBag;
P.Bag.OK     = zeros(1,P.Bag.n); % keeps track of nesting of Bags
P.Bag.iCavity= [];
P.Bag.iWall  = [];
P.Bag.iTube  = [];
P.Bag.iBag   = [];
for i=1:P.Bag.n
    BagIndex(i); % determines indices of bag-enclosed parts
end

% Baroreceptor location
P.Node.iBaro = Get('Node','Index',Net.Node.Baro);
% P.Patch.iPace= Get('Patch','Index',Net.Depolarization.Pace);
end


% ======== AUXILARY FUNCTIONS ============

function Element= DataWrite(ElementType,ElementName)
% 1st level
global P;
Element=P.(ElementType);
if isempty(Element.Name) % New P-structure
    Element.Name= ElementName;
    Element.n   = length(Element.Name);
else % existing P.Structure
    ReadIndex= Get(ElementType,'Index',ElementName);
    Element  = DataWrite1(Element,ReadIndex);
    Element.Name=Element.Name(ReadIndex);
    Element.n= length(ReadIndex);
end
end

function Element=DataWrite1(Element,ReadIndex) %ready for recursion
% removal of selected elements for all fields
Fields=fields(Element);
nF=length(Fields);
for iF=1:nF
    Field= Fields{iF};
    Aux  = Element.(Field);
    Incl = sum(strcmp(Field,ExcludedFields))==0;
    if isnumeric(Aux) && Incl && ~isempty(Aux)
        Element.(Field)=Element.(Field)(:,ReadIndex);
    end;
    if isstruct(Aux) % 1 step down for Adapt structure inside Element
        Element.(Field)=DataWrite1(Aux,ReadIndex);
    end
end
end

function str= ExcludedFields % excuded from DataWrite
str={'n','Adapt','ArtVen2NodeArt','ArtVen2NodeVen',...
    'Valve2NodeProx','Valve2NodeDist','AvValves','AvWalls',...
    'Tube2NodeProx','Tube2NodeDist',...
    'Cavity2Node','iBaro','iPace'...
    'c0','Valve2Node'}; %last 2 records should be removed in future
end


function BagIndex(iBag)
% function BagIndex(iBag)
% determines from indices of Chambers, TriSegs, ArtVens -> the indices of
% related Tubes, Walls and Cavities, enclosed by each Bag.

global P;
if P.Bag.OK(iBag)==1 % if indices already determined, skip this
    return; 
end;
iCh  = Get('Chamber','Index',P.Net.Bag.Chamber{iBag});
iTr  = Get('TriSeg' ,'Index',P.Net.Bag.TriSeg{iBag} );
iAv  = Get('ArtVen' ,'Index',P.Net.Bag.ArtVen{iBag} );
iCCh = P.Chamber.iCavity(iCh);
iCTr1= P.TriSeg.iCavity(iTr) ; iCTr=[iCTr1,iCTr1+1];
iCAv1= P.ArtVen.iCavity(iAv) ; iCAv=[iCAv1,iCAv1+1];
P.Bag.iCavity{iBag}= [iCCh,iCTr,iCAv];

iWCh = P.Chamber.iWall(iCh);
iWTr1= P.TriSeg.iWall(iTr) ; iWTr=[iWTr1,iWTr1+1,iWTr1+2];
P.Bag.iWall{iBag}=[iWCh,iWTr];

P.Bag.iTube{iBag}= NonZero(Get('Tube'   ,'Index',P.Net.Bag.Tube{iBag}   ));
jBag             = NonZero(Get('Bag'    ,'Index',P.Net.Bag.Bag{iBag}    ));
P.Bag.iBag{iBag} = jBag;
for j=1:length(jBag)
    jB=jBag(j);
    BagIndex(jB);
    P.Bag.iCavity{iBag}= [P.Bag.iCavity{iBag}, P.Bag.iCavity{jB}];
    P.Bag.iWall{iBag}  = [P.Bag.iWall{iBag}  , P.Bag.iWall{jB}  ];
    P.Bag.iTube{iBag}  = [P.Bag.iTube{iBag}  , P.Bag.iTube{jB}  ];
end
P.Bag.OK(iBag)=1;
end

function a=NonZero(a)
a=a(logical(a~=0));
end

