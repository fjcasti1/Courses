%% HOMEWORK 7 - FRANCISCO CASTILLO
%% Defined functions
%
% <include>./ADI_u.m</include>
%
% <include>./ADI_v.m</include>
%
% <include>./GaussTriSol.m</include>
%
% <include>./initialization.m</include>
%
% <include>./uvectors.m</include>
%
% <include>./vvectors.m</include>
%
% <include>./vvectors.m</include>
%
% <include>./poisson.m</include>
%
% <include>./multigrid.m</include>
%
% <include>./restrict.m</include>
%
% <include>./prolong.m</include>
%
% <include>./GaussSeidel.m</include>
%
% <include>./residual.m</include>
%
% <include>./InfNorm.m</include>
%
% <include>./divV.m</include>
%
% <include>./outletCorrection.m</include>
%
% <include>./LagProjection.m</include>
%

clear all; close all; format long; clc
axisSize=14;
markersize=16;
linewidth=3.5;
Lx = 4;
Ly = 2;
Mv = [1 0.5 0.25]*128;
Nv = [1 0.5 0.25]*64;
CFL = 1;
Re = 2;
nIterMax=50;
tol=1e-12;
p1(1)=0;
p2(1)=0;
p3(1)=0;
for i=1:length(Mv)
    M=Mv(i);
    N=Nv(i);
	if (M==128 && N==64)
        outputTime=[0.1 0.5 1 2];
    else
        outputTime=[0.1 0.5 1];
    end
    endtime=outputTime(end);
    hx = Lx/M;
    hy = Ly/N;
    if hx~=hy
        error('Cells not square')
    end
    time=0;
%     dt=CFL*0.25*hx^2*Re;
    dt=5e-3;
    n=1;
    % Define the points of the different meshes
    xu=linspace(0,4,M+1);
    yu=linspace(-hy/2,2+hy/2,N+2);
    xv=linspace(-hx/2,4+hx/2,M+2);
    yv=linspace(0,2,N+1);
    % Define initial values
    phi=zeros(M+2,N+2);
    [u,v]=initialization(M,N,hx,hy);
    iter=1;
    t=0;
    while time < endtime
        if (time < outputTime(n) && time+dt >= outputTime(n))
            dt=outputTime(n)-time;
            n=n+1;
        else
%             dt=CFL*0.25*hx^2*Re;
            dt=5e-3;
        end
        %% Solve for u(x,y,t)
        u=ADI_u(u,M,N,dt,hx,hy,Re);
        %% Solve for v(x,y,t)
        v=ADI_v(v,M,N,dt,hx,hy,Re);
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
        %% Filled contour plots, only for M=128 and N=64
        if (M==128 && N==64 && ismember(time,outputTime))
            % Plot u
            figure(n-1)
            contourf(xu,yu,u')
            hold on
            plot(1,0.5,'r.','markersize',markersize,'linewidth',linewidth)
            axis([0 4 0 2])
            caxis([-1.5 3])
            colorbar
            xlabel('$x$','Interpreter','latex')
            ylabel('$y$','Interpreter','latex')
            yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
            xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
            set(get(gca,'ylabel'),'rotation',0)
            set(gca,'fontsize',axisSize)
            pbaspect([2 1 1])
            grid on
            txt=['Latex/FIGURES/u_' num2str(n-1)];
            saveas(gcf,txt,'epsc')
            % Plot v
            figure(n+3)
            contourf(xv,yv,v')
            hold on
            plot(1,1.5,'r.','markersize',markersize,'linewidth',linewidth)
            axis([0 4 0 2])
            caxis([-3 1])
            colorbar
            xlabel('$x$','Interpreter','latex')
            ylabel('$y$','Interpreter','latex')
            yticks([0 0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0]);
            xticks([0 0.50 1.0 1.50 2.0 2.50 3.0 3.50 4.0]);
            set(get(gca,'ylabel'),'rotation',0)
            set(gca,'fontsize',axisSize)
            pbaspect([2 1 1])
            grid on
            txt=['Latex/FIGURES/v_' num2str(n-1)];
            saveas(gcf,txt,'epsc')
        end
        %% Probes, only for M=128 and N=64
        if (M==128 && N==64 && time<=2+dt)
            iter=iter+1;
            t=[t;time];
            p1(iter) = (u (find(xu==1),find(yu<=0.5,1,'last'))...
                +u (find(xu==1),find(yu>=0.5,1)))/2;
            p2(iter) = (v (find(xv<=1,1,'last'),find(yv==1.5))...
                +v (find(xv>=1,1),find(yv==1.5)))/2;
        end
        %% Probes value for u, v and Y at t=1
        if time==1
            p11(i)= (u (find(xu==1),find(yu<=0.5,1,'last'))...
                    +u (find(xu==1),find(yu>=0.5,1)))/2;
            p21(i)= (v (find(xv<=1,1,'last'),find(yv==1.5))...
                    +v (find(xv>=1,1),find(yv==1.5)))/2;
        end
    end
    %% Plot the probes, only for M=128 and N=64
    if (M==128 && N==64)
        figure
        plot(t,p1,'linewidth',2)
        xlim([0 2])
        xlabel('$t$','Interpreter','latex')
        ylabel('$u(1,0.5,t)$','Interpreter','latex')
        set(gca,'fontsize',axisSize)
        grid on
        txt='Latex/FIGURES/probe1';
        saveas(gcf,txt,'epsc')
        figure
        plot(t,p2,'linewidth',2)
        xlim([0 2])
        xlabel('$t$','Interpreter','latex')
        ylabel('$v(1,1.5,t)$','Interpreter','latex')
        set(gca,'fontsize',axisSize)
        grid on
        txt='Latex/FIGURES/probe2';
        saveas(gcf,txt,'epsc')
    end
end
%% GCI analysis 
clc
r=2
Fsec=1.25
% Probe 1
p_u=log((p11(3)-p11(2))/(p11(2)-p11(1)))/log(r)
u_h0=p11(1)+(p11(1)-p11(2))/(r^p_u-1)
GCI21_u=Fsec*(p11(2)-p11(1))/(p11(1)*(r^p_u-1))
GCI32_u=Fsec*(p11(3)-p11(2))/(p11(2)*(r^p_u-1))
coeff_u=GCI21_u*r^p_u/GCI32_u
percent_u=GCI21_u*100
pause

% Probe 2
p_v=log((p21(3)-p21(2))/(p21(2)-p21(1)))/log(r)
v_h0=p21(1)+(p21(1)-p21(2))/(r^p_v-1)
GCI21_v=Fsec*(p21(2)-p21(1))/(p21(1)*(r^p_v-1))
GCI32_v=Fsec*(p21(3)-p21(2))/(p21(2)*(r^p_v-1))
coeff_v=GCI21_v*r^p_v/GCI32_v
percent_v=GCI21_v*100
% pause
