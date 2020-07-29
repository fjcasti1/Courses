function Problem2(M,N,CFL,boolMovie,design)
%% FINAL PROJECT - FRANCISCO CASTILLO
%% Defined functions
%
%


close all; format long; clc
% delete(gcp('nocreate'))
% poolWorkers = parpool('local',2);
tic
%% Create directories to store results
txt=['RESULTS_CFL0',num2str(10*CFL)];
path=[txt,'/Results_M',num2str(M)];
pathGCI=[txt,'/GCIdata_M',num2str(M)];
pathFIGURES=[txt,'/Results_M',num2str(M),'/FIGURES/'];
pathMOVIES=[txt,'/Results_M',num2str(M),'/MOVIES/'];
mkdir(pathFIGURES)
if boolMovie
	mkdir(pathMOVIES)
end
%% Global Variables
global Lx;
global Ly;
global hx; global hy;
global xu; global yu;
global xv; global yv;
global xY; global yY;
global Re; global Sc;
global a1; global b1; global uavg1; global startTime1;
global a2; global b2; global uavg2; global startTime2;
global a3; global b3; global uavg3; global startTime3;
global ao; global bo;
global nIterMax; global tol;

%% Read parameters
DATA=importdata(['Designs/Design',num2str(design),'.txt'],',');
designData=DATA.data;
a1=designData(1,1);
b1=designData(1,2);
uavg1=designData(1,3);
startTime1=designData(1,4);
a2=designData(2,1);
b2=designData(2,2);
uavg2=designData(2,3);
startTime2=designData(2,4);
a3=designData(3,1);
b3=designData(3,2);
uavg3=designData(3,3);
startTime3=designData(3,4);
ao=designData(5,1);
bo=designData(5,2);

%% Mesh parameters
Lx=2;
Ly=4;
hx = Lx/M;
hy = Ly/N;

xu=linspace(0,Lx,M+1);
yu=linspace(-hy/2,Ly+hy/2,N+2);
xv=linspace(-hx/2,Lx+hx/2,M+2);
yv=linspace(0,Ly,N+1);
xY=linspace(-5*hx/2,Lx+5*hx/2,M+6);
yY=linspace(-5*hx/2,Ly+5*hx/2,N+6);
%% Other parameters
Re = 50;
Sc=1;
nIterMax=20;
tol=1e-5;

outputTime=[1 3 5 10 15 20];
endtime=outputTime(end);
n=1;
tstep=1;
movieTime=0.05:0.05:20;
nMovie=1;

if hx~=hy
    error('Cells not square')
end

disp(' ')
if boolMovie
disp(['--> STARTED Movie Simulation with M=',num2str(M),', N=',num2str(N),', CFL=',num2str(CFL),' and design number ',num2str(design)])
else
disp(['--> STARTED Simulation with M=',num2str(M),', N=',num2str(N),', CFL=',num2str(CFL),' and design number ',num2str(design)])
end
disp(' ')
pause(5)

time=0;
% Define initial values

[u,v,Y]=initialization(M,N);

if boolMovie
	getMovieFrame(u,v,Y,nMovie,pathMOVIES,time)
end

phi=zeros(M+2,N+2);
dt=stableTimeStep(u,v,CFL);

Hu=zeros(size(u));
Hv=zeros(size(v));
HY=zeros(size(Y));

R=zeros(5e2,1);
S=zeros(5e2,1);
t=zeros(5e2,1);

[R(tstep,1),S(tstep,1)]=performance(Y,v,M,N);

while time < endtime
    %% Calculate Time Step
    if boolMovie
        [dt,nMovie]=calcTimeStep(time,dt,movieTime(nMovie),nMovie,u,v,CFL);
    else
        dt=calcTimeStep(time,dt,outputTime(n),n,u,v,CFL);
    end
    %% Solve Hyperbolic part of Y equation
    HY=HyperbolicY(Y,M,N,time,dt,u,v,HY);
    %% Solve Y equation 
    Y=ADI_Y(Y,M,N,dt,time,HY);
    %% Solve Burger's viscous equation
    [u,v,Hu,Hv]=solveBurgers2D(u,v,Hu,Hv,M,N,dt,time);
    %% Outlet correction
    v=outletCorrection(u,v);
    %% Calculate rhs
    rhs=divV(u,v,M,N,dt);
    %% Solve Poisson equation, V-cycle multigrid
    phi = poisson(phi,rhs);
    %% Project velocities using Lagrange multiplier
    [u,v]=LagProjection(u,v,M,N,phi,dt);
    %% Next time step
    time=time+dt;
    tstep=tstep+1;
    %% Store performance parameters
    if mod(tstep,5e2)==0
        R=[R;zeros(5e2,1)];
        S=[S;zeros(5e2,1)];
        t=[t;zeros(5e2,1)];
        elapsedTime=seconds(toc);
        [hours,mins,secs]=hms(elapsedTime);
        if boolMovie
            disp(['Movie Simulation with M=',num2str(M),', N=',num2str(N),', CFL=',num2str(CFL),' and design number ',num2str(design),...
            ' has performed ',num2str(tstep),' iterations. Time step=',num2str(dt),'. Elapsed time is ',num2str(hours),' hours, ',num2str(mins),...
            ' minutes and ',num2str(secs),' seconds.'])
        else
            disp(['Simulation with M=',num2str(M),', N=',num2str(N),', CFL=',num2str(CFL),' and design number ',num2str(design),...
            ' has performed ',num2str(tstep),' iterations. Time step=',num2str(dt),'. Elapsed time is ',num2str(hours),' hours, ',num2str(mins),...
            ' minutes and ',num2str(secs),' seconds.'])
        end
    end
    [R(tstep,1),S(tstep,1)]=performance(Y,v,M,N);
    t(tstep)=time; % time vector used to plot

    %% Save progress
    if ismember(time,outputTime)
        % Check divergence free
        rhscheck=InfNorm(divV(u,v,M,N,dt));
        elapsedTime=seconds(toc);
        [hours,mins,secs]=hms(elapsedTime);
        if boolMovie
            disp(['Movie Simulation with M=',num2str(M),', N=',num2str(N),', CFL=',num2str(CFL),' and design number ',num2str(design),...
            ' has reached time=',num2str(time),', and having a InfNorm(rhs)=',num2str(rhscheck),'. Elapsed time is ',num2str(hours),' hours, ',num2str(mins),...
            ' minutes and ',num2str(secs),' seconds.'])
        else
            disp(['Simulation with M=',num2str(M),', N=',num2str(N),', CFL=',num2str(CFL),' and design number ',num2str(design),...
            ' has reached time=',num2str(time),', and having a InfNorm(rhs)=',num2str(rhscheck),'. Elapsed time is ',num2str(hours),' hours, ',num2str(mins),...
            ' minutes and ',num2str(secs),' seconds.'])
        end
        txt=[path,'_t',num2str(time)];
        save(txt)
    end
    %% Frames capture
    if ( boolMovie && ismember(time,movieTime) )
        getMovieFrame(u,v,Y,nMovie,pathMOVIES,time)
    end
    %% Output Contour Plots
    if ismember(time,outputTime)
        ContourPlots(u,v,Y,n,pathFIGURES)
        n=n+1;
    end
close all
end
%% Output Curve Plots
CurvePlots(R,S,t,pathFIGURES)
%% Calculate T
T=calculateT(S,t);
%% Save GCI data
elapsedTime=seconds(toc);
[hours,mins,secs]=hms(elapsedTime);
save(pathGCI,'M','N','CFL','T','hours','mins','secs')
disp(' ')
if boolMovie
    disp(['--> FINISHED Movie Simulation with M=',num2str(M),', N=',num2str(N),', CFL=',num2str(CFL),' and design number ',num2str(design),...
            ' obtaining T=',num2str(T)])
else
    disp(['--> FINISHED Simulation with M=',num2str(M),', N=',num2str(N),', CFL=',num2str(CFL),' and design number ',num2str(design),...
            ' obtaining T=',num2str(T)])
end
disp(['Elapsed time is ',num2str(hours),' hours, ',num2str(mins),' minutes and ',num2str(secs),' seconds.'])
disp(' ')

end
