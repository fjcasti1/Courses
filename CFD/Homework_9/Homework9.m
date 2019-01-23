%% HOMEWORK 9 - FRANCISCO CASTILLO
%% Defined functions
%
% <include>./EntropyFix.m</include>
%
% <include>./HY_alpha.m</include>
%
% <include>./HY_beta.m</include>
%
% <include>./HY_G.m</include>
%
% <include>./HY_Phi.m</include>
%
% <include>./initialCondition.m</include>
%
% <include>./numericalFlux.m</include>
%
% <include>./sigmaG.m</include>
%
% <include>./TVD2order.m</include>
%
%% ProblemHomework8.m
clear variables
close all
clc
format long

axisSize=14;
linewidth=1.5;
L=4;
CFL=0.8;
M=32;
i=0;
T=nan(3,7);
check=1;
while check>0.01
    i=i+1;
    M=2*M
    dx=L/M;
    x=linspace(3-1.5*dx,7+1.5*dx,M+4)'; % Cell centered mesh with two
                                        % ghost cells at each side                         
    %% Initialization
    time=0;
    u=initialCondition(x,M);
    % Plot initial condition
    if (M==256 || M==1024)
            figure(1)
            plot(x,u,'linewidth',linewidth)
            grid on
            axis([min(x) max(x) -1 3])
            xlabel('$x$','Interpreter','latex')
            ylabel('$u(x,t)$','Interpreter','latex')
            set(gca,'fontsize',axisSize)
            if M==256
                txt='Latex/FIGURES/uinitial_M256';
            elseif M==1024
                txt='Latex/FIGURES/uinitial_M1024';
            end
            saveas(gcf,txt,'epsc')
    end
    E=0.5*u.^2;
    alpha=HY_alpha(u,E,M);
    dt=CFL*dx/(max(abs(alpha)));
    
	outputTime=[0.1 0.5 1 3.0];
    endtime=outputTime(end);
    
    n=1;
    while time < endtime
        E=0.5*u.^2;
        alpha=HY_alpha(u,E,M);
        if (time < outputTime(n) && time+dt >= outputTime(n))
            dt=outputTime(n)-time;
            n=n+1;
        else
            dt=CFL*dx/(max(abs(alpha)));
        end
        G=HY_G(u,alpha,dt,dx,M);
        beta=HY_beta(u,G,M);
        psi=EntropyFix(alpha(2:end-1)+beta);
        Phi=HY_Phi(u,G,alpha,beta,M);
        h = numericalFlux(E,Phi,M);
        
        u=TVD2order(u,h,dt,dx,M);
        time=time+dt;
        %% Plot required at part 5
        if (M==256 && ismember(time,outputTime))
            figure(n)
            plot(x,u,'linewidth',linewidth)
            grid on
            axis([3 7 -1 3])
            xlabel('$x$','Interpreter','latex')
            ylabel('$u(x,t)$','Interpreter','latex')
            set(gca,'fontsize',axisSize)
            txt=['Latex/FIGURES/u_M256_' num2str(n-1)];
            saveas(gcf,txt,'epsc')
        end
        if (M==1024 && ismember(time,outputTime))
            figure(n)
            plot(x,u,'linewidth',linewidth)
            grid on
            axis([3 7 -1 3])
            xlabel('$x$','Interpreter','latex')
            ylabel('$u(x,t)$','Interpreter','latex')
            set(gca,'fontsize',axisSize)
            txt=['Latex/FIGURES/u_M1024_' num2str(n-1)];
            saveas(gcf,txt,'epsc')
        end
        %% u at x=6 and t=1
        if time==1
            uGCI(i)=(u(find(x<=6,1,'last'))+u(find(x>=6,1)))/2;
        end
    end
    %% GCI analysis 
    if i>=3
        r=2;
        Fsec=1.25;
        p(i)=log(abs(uGCI(i-2)-uGCI(i-1))/abs(uGCI(i-1)-uGCI(i)))/log(r);
        u_h0(i)=uGCI(i)+(uGCI(i)-uGCI(i-1))/(r^p(i)-1);
        GCI12(i)=Fsec*abs(1-uGCI(i-1)/uGCI(i))/(r^p(i)-1);
        GCI23(i)=Fsec*abs(1-uGCI(i-2)/uGCI(i-1))/(r^p(i)-1);
        coeff(i)=GCI12(i)*r^p(i)/GCI23(i);
        percent(i)=GCI12(i)*100;
        check=abs(percent(i));
        % Include results in a table
        T(i,2)=p(i);
        T(i,3)=u_h0(i);
        T(i,4)=GCI12(i);
        T(i,5)=GCI23(i);
        T(i,6)=coeff(i);
        T(i,7)=percent(i);
    end
    % Include results in a table
    T(i,1)=M;
end
T=array2table(T,'VariableNames',{'M','p','phi0','GCI12','GCI23','coeff','Check'})
if CFL==0.1
    save('Case1_CFL01');
elseif CFL==0.5
    save('Case2_CFL05');
end