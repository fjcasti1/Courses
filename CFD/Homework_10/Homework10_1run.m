function Homework10_1run(M,N)
%% HOMEWORK 10 - FRANCISCO CASTILLO

%% Defined functions
%
%
close all; format long; clc

Lx=4;
Ly=2;
% M=64;
% N=32;
CFL = 0.8;
Re = 20;
Sc=0.5;

% put5=zeros(6,1);
% pvt5=zeros(6,1);
% pYt5=zeros(6,1);
% Rt5=zeros(6,1);
% 
% put10=zeros(6,1);
% pvt10=zeros(6,1);
% pYt10=zeros(6,1);
% Rt10=zeros(6,1);

% i=0;
% checku5=10;
% checkv5=10;
% checkY5=10;

outputTime=[1 2 3 5 7 10 13 16 19 22 25 28 30];
endtime=outputTime(end);
    
% while ( checku5>0.02 || checkv5>0.02 || checkY5>0.4 )
% while ( checku>1 && checkv>1.2 && checkY>6 )
%     i=i+1;
%     i=1;
%     M=2*M
%     N=2*N;
    hx = Lx/M;
    hy = Ly/N;
    if hx~=hy
        error('Cells not square')
    end
    time=0
    n=1;
    tstep=1;
    
    % Define initial values
    [u,v,Y]=initialization(M,N,hx,hy);
    dt=calcTimeStep(u,v,hx,hy,CFL);
    
    Hu=zeros(size(u));
    Hv=zeros(size(v));
    
    pu=zeros(5e2,1);
    pv=zeros(5e2,1);
    pY=zeros(5e2,1);
    R=zeros(5e2,1);
    t=zeros(5e2,1);
    
	[pu(tstep,1),pv(tstep,1),pY(tstep,1)]=probeValues(u,v,Y,M,N,hx,hy);
    R(tstep,1)=performance(Y,M,N,hx,hy,Lx,Ly);
    
    while time < endtime
        %% Calculate Time Step
        if (time < outputTime(n) && time+dt >= outputTime(n))
            dt=outputTime(n)-time;
            n=n+1;
        else
            dt=calcTimeStep(u,v,hx,hy,CFL);
        end
        %% Solve Hyperbolic part of Y equation
        HY=HyperbolicY(Y,M,N,hx,hy,dt,u,v);
        %% Solve Y equation 
        Y=ADI_Y(Y,M,N,dt,hx,hy,Re,Sc,HY);
        %% Solve Burger's viscous equation
        [u,v,Hu,Hv]=solveBurgers2D(u,v,Hu,Hv,M,N,hx,hy,dt,Re,time);
        %% Next time step
        time=time+dt;
        tstep=tstep+1;
        %% Store probe values and performance parameter
        if mod(tstep,5e2)==0
            pu=[pu;zeros(5e2,1)];
            pv=[pv;zeros(5e2,1)];
            pY=[pY;zeros(5e2,1)];
            R=[R;zeros(5e2,1)];
            t=[t;zeros(5e2,1)];
        end
        [pu(tstep,1),pv(tstep,1),pY(tstep,1)]=probeValues(u,v,Y,M,N,hx,hy);
        R(tstep,1)=performance(Y,M,N,hx,hy,Lx,Ly);
        t(tstep)=time; % time vector used to plot
        if time==5
            put5=pu(tstep);
            pvt5=pv(tstep);
            pYt5=pY(tstep);
            Rt5=R(tstep);
        end
        if time==endtime
            putend=pu(tstep);
            pvtend=pv(tstep);
            pYtend=pY(tstep);
            Rtend=R(tstep);
        end
%         %% Outputs
%         if ( M==256 && N==128 )
%             if ismember(time,outputTime)
%                 ContourPlots(u,v,Y,M,N,hx,hy,n-1)
%             end
%             if time==endtime
%                 CurvePlots(pu,pv,pY,R,t)
%             end
%         end
        %% Save progress
        if ismember(time,outputTime)
            time
            txt=['DataSavedSteady/Save_M',num2str(M),'_time',num2str(time)];
            save(txt)
        end
    end
    
    %% Save GCI data
    txt=['DataSavedSteady/GCIdata_Steady/GCIdata_M',num2str(M)];
	save(txt,'M','N','put5','pvt5','pYt5','Rt5','putend','pvtend','pYtend','Rtend')
end