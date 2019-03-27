function OutNew=SteadyStateP
% function OutNew=SteadyStateP
% OutNew= new estimate vector of steady state parameters, to be used 
% as start condition for next heart beat
% Theo Arts, Maastricht University, Oct 13, 2012

global P

% read input-output record, take logarithm for relative change
In  = log(P.General.In );
Out = log(P.General.Out);
[nx,np]=size(In); % [number of samples, number of parameters]
yRef= Out(end,:);

nxn= 1+mod(np,5);%ceil(2+mod(np,4));
nMin= max(2,nx-nxn); % number of samples used for prediction
Rgn=(nMin:nx)'; %indices used for prediction
nxn=length(Rgn);

if nxn>1
    X=In(Rgn,:); Y=Out(Rgn,:); %U=Out(Rgn-1,:);
    nt=size(X,1); Col1=ones(nt,1);
    YmX =Y-X;
    W1 = 1.0 ./ (0.01+sqrt(mean(YmX.^2,2)));
    %     W1=0*W1+1;
    W=repmat(W1,[1,np]);
    M1 = [W.*YmX ,W1.*Col1];
    N1 = pinv(M1) * (W.*(Y-Col1*yRef));
    dy = N1(end,:);
    y1=yRef+dy; % reasonable estimate
    
    %Limit Eigen values of dy/dx < 1.
    DX=X-repmat(y1,[nt,1]);
    DY=Y-repmat(y1,[nt,1]);
    M=pinv(DX)*DY;
    [p1 q1 r1]=svd(M); a=1.1;
    q1=a*tanh(q1/a);%min(q1,1.0); % Eigen values <1.0
    M2=p1*q1*r1';
    DX=X-repmat(yRef,[nt,1]);
    DY=Y-repmat(yRef,[nt,1]);
    dy=mean((DY-DX*M2)*pinv(eye(np)-M2));
    yNew=yRef+dy; % more stable estimate
else
    yNew=yRef;
end
    
OutNew=exp(yNew); % prediction after exponential

end

