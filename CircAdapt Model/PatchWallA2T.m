function PatchWallA2T
% function PatchWallA2T
% Determines linear function patch area Am = fu Tension T
% Am(T)=Am0+T*DADT
% Constants Am0 and DADT are determined from Patch state variables
% and reference geometry
% Theo Arts, Maastricht University, Oct 13, 2012

global P;

% Patch Am= Am0+DADT*T, provides Am0 and DADT
AmRef= P.Patch.AmRef; % midwall area for ref sarcomere length 2mu
LsRef= P.Patch.LsRef; % ref sarcomere length 2mu
VWall= P.Patch.VWall; % wall volume

Lsi   = P.Patch.Lsi; % unloaded sarc length= state variable
Lambda= bsxfun(@rdivide,Lsi      ,LsRef); %extension factor
Am    = bsxfun(@times  ,Lambda.^2,AmRef); %actual midwall area
Ef    = log(Lambda); %natural fiber strain

P.Patch.Ef= Ef;
SarcEf2Sf; % sarcomere strain->stress
Sf   = P.Patch.Sf; % sarcomere stress
T    = bsxfun(@times  , Sf, (0.5 *VWall)) ./ Am   ; % tension
DADT = bsxfun(@rdivide, Am.^2 ./ max(P.Patch.DSfDEf-2*Sf,1e3),...
    0.25*VWall); %compliance
P.Patch.T   = T;
P.Patch.DADT= DADT;
P.Patch.Am0 = Am - T.*DADT; % zero load midwall area

% Wall is composed of patches: Also for wall: Am(T)=Am0+DADT*T
for iWall=1:P.Wall.n
    iPatch= (P.Wall.iPatch(iWall)-1)+(1:P.Wall.nPatch(iWall));
    P.Wall.VWall(iWall)  = sum(P.Patch.VWall(iPatch));
    P.Wall.Am0(:,iWall)  = P.Wall.AmDead(iWall)+sum(P.Patch.Am0(:,iPatch),2);
    P.Wall.DADT(:,iWall) = sum(P.Patch.DADT(:,iPatch),2);
end
end

