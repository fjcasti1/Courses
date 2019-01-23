%% HOMEWORK 10 - FRANCISCO CASTILLO
%% Defined functions
%
% <include>./initialization.m</include>
%
% <include>./applyBCs.m</include>
%
% <include>./calcTimeStep.m</include>
%
% <include>./ADI_u.m</include>
%
% <include>./uvectors.m</include>
%
% <include>./ADI_v.m</include>
%
% <include>./vvectors.m</include>
%
% <include>./ADI_Y.m</include>
%
% <include>./Yvectors.m</include>
%
% <include>./HyperbolicY.m</include>
%
% <include>./TVDRK3_2D.m</include>
%
% <include>./WENO5_2D.m</include>
%
% <include>./psiWENO.m</include>
%
% <include>./HyperbolicBurgers2D.m</include>
%
% <include>./solveBurgers2D.m</include>
%
% <include>./performance.m</include>
%
% <include>./GaussTriSol.m</include>
%
% <include>./ContourPlots.m</include>
%
% <include>./probeValues.m</include>
%
% <include>./CurvePlots.m</include>
%
%%
clear all; close all; format long; clc

Lx=4;
Ly=2;
M=64;
N=32;
CFL = 0.8;
Re = 20;
Sc=0.5;

put5=zeros(6,1);
pvt5=zeros(6,1);
pYt5=zeros(6,1);
Rt5=zeros(6,1);

put10=zeros(6,1);
pvt10=zeros(6,1);
pYt10=zeros(6,1);
Rt10=zeros(6,1);

Tu5=nan(5,7);
Tv5=nan(5,7);
TY5=nan(5,7);

Tu10=nan(5,7);
Tv10=nan(5,7);
TY10=nan(5,7);


i=0;
checku5=10;
checkv5=10;
checkY5=10;

outputTime=[1 2 3 5 7 10];
endtime=outputTime(end);
    
while ( checku5>0.02 || checkv5>0.02 || checkY5>0.4 )
% while ( checku>1 && checkv>1.2 && checkY>6 )
    i=i+1;
    M=2*M
    N=2*N;
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
            put5(i)=pu(tstep);
            pvt5(i)=pv(tstep);
            pYt5(i)=pY(tstep);
            Rt5(i)=R(tstep);
        end
        if time==10
            put10(i)=pu(tstep);
            pvt10(i)=pv(tstep);
            pYt10(i)=pY(tstep);
            Rt10(i)=R(tstep);
        end
        %% Outputs
        if ( M==256 && N==128 )
            if ismember(time,outputTime)
                ContourPlots(u,v,Y,M,N,hx,hy,n-1)
            end
            if time==endtime
                CurvePlots(pu,pv,pY,R,t)
            end
        end
        %% Save progress
        if ismember(time,outputTime)
            time
            txt=['DataSaved/Save_M',num2str(M),'_time',num2str(time)];
            save(txt)
        end
    end
    %% GCI analysis 
    % Values at t=5
    Tu5(i,1)=M;  Tu5(i,2)=N;
    Tv5(i,1)=M;  Tv5(i,2)=N;
    TY5(i,1)=M;  TY5(i,2)=N;
    if i>=3
        [Tu5,Tv5,TY5]=GCIanalysis(put5,pvt5,pYt5,i,Tu5,Tv5,TY5);
        checku5=abs(Tu5(i,5));
        checkv5=abs(Tv5(i,5));
        checkY5=abs(TY5(i,5));
    end
    % Values at t=10
    Tu10(i,1)=M;  Tu10(i,2)=N;
    Tv10(i,1)=M;  Tv10(i,2)=N;
    TY10(i,1)=M;  TY10(i,2)=N;
    if i>=3
        [Tu10,Tv10,TY10]=GCIanalysis(put10,pvt10,pYt10,i,Tu10,Tv10,TY10);
        checku10=abs(Tu10(i,5));
        checkv10=abs(Tv10(i,5));
        checkY10=abs(TY10(i,5));
    end
end

Tu5=array2table(Tu5,'VariableNames',{'M','N','p','phi0','GCI12','GCI23','beta'})
Tv5=array2table(Tv5,'VariableNames',{'M','N','p','phi0','GCI12','GCI23','beta'})
TY5=array2table(TY5,'VariableNames',{'M','N','p','phi0','GCI12','GCI23','beta'})

Tu10=array2table(Tu10,'VariableNames',{'M','N','p','phi0','GCI12','GCI23','beta'})
Tv10=array2table(Tv10,'VariableNames',{'M','N','p','phi0','GCI12','GCI23','beta'})
TY10=array2table(TY10,'VariableNames',{'M','N','p','phi0','GCI12','GCI23','beta'})