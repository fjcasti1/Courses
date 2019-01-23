function Problem1(M,N,CFL,boolMovie)
%% FINAL PROJECT - FRANCISCO CASTILLO
%% Defined functions
%
%

close all; format long; clc
if nargin<4
    boolMovie=0;
end
tic
%% Create directories to store results
txt=['RESULTS_CFL0',num2str(10*CFL)];
path=[txt,'/Results_M',num2str(M),'/Results_M',num2str(M)];
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
global hx;
global hy;
global xu;
global yu;
global xv;
global yv;
global xY;
global yY;
global Re;
global Sc;

Lx=4;
Ly=2;
hx = Lx/M;
hy = Ly/N;

xu=linspace(0,4,M+1);
yu=linspace(-hy/2,2+hy/2,N+2);
xv=linspace(-hx/2,4+hx/2,M+2);
yv=linspace(0,2,N+1);
xY=linspace(-5*hx/2,4+5*hx/2,M+6);
yY=linspace(-5*hx/2,2+5*hx/2,N+6);
% CFL = 0.8;
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
disp(['--> STARTED Movie Simulation with M=',num2str(M),' and N=',num2str(N),' and CFL=',num2str(CFL)])
else
disp(['--> STARTED Simulation with M=',num2str(M),' and N=',num2str(N),' and CFL=',num2str(CFL)])
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
    HY=HyperbolicY(Y,M,N,dt,u,v,HY);
    %% Solve Y equation 
    Y=ADI_Y(Y,M,N,dt,HY);
    %% Solve Burger's viscous equation
    [u,v,Hu,Hv]=solveBurgers2D(u,v,Hu,Hv,M,N,dt,time);
    %% Outlet correction
    v=outletCorrection(u,v,xv,hx,hy);
    %% Calculate rhs
    rhs=divV(u,v,M,N,hx,hy,dt);
    %% Solve Poisson equation, V-cycle multigrid
    phi = poisson(phi,rhs,hx,hy,nIterMax,tol);
    %% Project velocities using Lagrange multiplier
    [u,v]=LagProjection(u,v,phi,dt,hx,hy);
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
            disp(['Movie Simulation with M=',num2str(M),' and N=',num2str(N),' and CFL=',num2str(CFL)...
            ' has performed ',num2str(tstep),' iterations. Time step=',num2str(dt)])
            disp(['--> Elapsed time is ',num2str(hours),' hours, ',num2str(mins),' minutes and ',num2str(secs),' seconds.'])
        else
            disp(['Simulation with M=',num2str(M),' and N=',num2str(N),' and CFL=',num2str(CFL),...
            ' has performed ',num2str(tstep),' iterations. Time step=',num2str(dt)])
            disp(['--> Elapsed time is ',num2str(hours),' hours, ',num2str(mins),' minutes and ',num2str(secs),' seconds.'])
        end
    end
    [R(tstep,1),S(tstep,1)]=performance(Y,v,M,N);
    t(tstep)=time; % time vector used to plot

    %% Save progress
    if ismember(time,outputTime)
        % Check divergence free
        rhscheck=InfNorm(divV(u,v,M,N,hx,hy,dt));
        if boolMovie
            disp(['Movie Simulation with M=',num2str(M),' and N=',num2str(N),' and CFL=',num2str(CFL),...
            ' has reached time=',num2str(time),', and having a InfNorm(rhs)=',num2str(rhscheck)])
        else
            disp(['Simulation with M=',num2str(M),' and N=',num2str(N),' and CFL=',num2str(CFL),...
            ' has reached time=',num2str(time),', and having a InfNorm(rhs)=',num2str(rhscheck)])
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
save(pathGCI,'M','N','CFL','T')
disp(' ')
if boolMovie
    disp(['--> FINISHED Movie Simulation with M=',num2str(M),' and N=',num2str(N),' and CFL=',num2str(CFL),...
            ' obtaining T=',num2str(T)])
else
    disp(['--> FINISHED Simulation with M=',num2str(M),' and N=',num2str(N),' and CFL=',num2str(CFL),...
            ' obtaining T=',num2str(T)])
end
elapsedTime=seconds(toc);
[hours,mins,secs]=hms(elapsedTime);
disp(['--> Elapsed time is ',num2str(hours),' hours, ',num2str(mins),' minutes and ',num2str(secs),' seconds.'])
disp(' ')

end
