%% HOMEWORK 8 - FRANCISCO CASTILLO
%% Defined functions
%
% <include>./velocity.m</include>
%
% <include>./updateGhostCells.m</include>
%
% <include>./phi_left.m</include>
%
% <include>./WENO5.m</include>
%
% <include>./psiWENO.m</include>
%
% <include>./TVDRK3.m</include>
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
M=64;
i=0;
T=nan(3,7);
check=1;
while check>0.1
    i=i+1;
    M=2*M
    hx=L/M;
    x=linspace(-1-2.5*hx,3+2.5*hx,M+6)'; % Cell centered mesh with three
                                        % ghost cells at each side
    %% Initialization
    time=0;
    a=velocity(x,time);
    phi = zeros(M+6,1);
    phi(3)=2*phi_left(time); 
    phi(2)=2*phi_left(time);  % Since phi(4:6)=0 it is not necessary
    phi(1)=2*phi_left(time);  % to include them

    dt=CFL*hx/(max(abs(a)));
    outputTime=[0.25 0.5 1 1.25 1.5 2.1];
    endtime=outputTime(end);

    n=1;
    step=0;
    while time < endtime
        a=velocity(x,time);
        if (time < outputTime(n) && time+dt >= outputTime(n))
            dt=outputTime(n)-time;
            n=n+1;
        else
            dt=CFL*hx/(max(abs(a)));
        end
        phi = TVDRK3(phi,M,hx,dt,time,a);
        time=time+dt;
        %% Plot required at part 5
        if (M==256 && CFL==0.8 && ismember(time,outputTime))
            figure(n-1)
            plot(x,phi,'linewidth',linewidth)
            grid on
            axis([-1 3 -0.1 1.1])
            xlabel('$x$','Interpreter','latex')
            ylabel('$\phi(x,t)$','Interpreter','latex')
            set(gca,'fontsize',axisSize)
            txt=['Latex/FIGURES/phi_' num2str(n-1)];
            saveas(gcf,txt,'epsc')
        end
        %% Plot required at part 6
        % For CFL=0.8
        if (M==4096 && CFL==0.8 && ismember(time,outputTime))
            figure(n+6)
            plot(x,phi,'linewidth',linewidth)
            grid on
            axis([-1 3 -0.1 1.1])
            xlabel('$x$','Interpreter','latex')
            ylabel('$\phi(x,t)$','Interpreter','latex')
            set(gca,'fontsize',axisSize)
            txt=['Latex/FIGURES/phi08_' num2str(n-1)];
            saveas(gcf,txt,'epsc')
        end
        % For CFL=0.5
        if (M==1024 && CFL==0.5 && ismember(time,outputTime))
            figure(n+6)
            plot(x,phi,'linewidth',linewidth)
            grid on
            axis([-1 3 -0.1 1.1])
            xlabel('$x$','Interpreter','latex')
            ylabel('$\phi(x,t)$','Interpreter','latex')
            set(gca,'fontsize',axisSize)
            txt=['Latex/FIGURES/phi05_' num2str(n-1)];
            saveas(gcf,txt,'epsc')
        end
        %% Phi at x=0 and t=1.25
        if time==1.25
            phiGCI(i)=(phi(find(x<=0,1,'last'))+phi(find(x>=0,1)))/2;
        end
    end
    %% GCI analysis 
    if i>=3
        r=2;
        Fsec=1.25;
        p(i)=log(abs(phiGCI(i-2)-phiGCI(i-1))/abs(phiGCI(i-1)-phiGCI(i)))/log(r);
        phi_h0(i)=phiGCI(i)+(phiGCI(i)-phiGCI(i-1))/(r^p(i)-1);
        GCI12(i)=Fsec*abs(1-phiGCI(i-1)/phiGCI(i))/(r^p(i)-1);
        GCI23(i)=Fsec*abs(1-phiGCI(i-2)/phiGCI(i-1))/(r^p(i)-1);
        coeff(i)=GCI12(i)*r^p(i)/GCI23(i);
        percent(i)=GCI12(i)*100;
        check=abs(percent(i));
        % Include results in a table
        T(i,2)=p(i);
        T(i,3)=phi_h0(i);
        T(i,4)=GCI12(i);
        T(i,5)=GCI23(i);
        T(i,6)=coeff(i);
        T(i,7)=percent(i);
    end
    % Include results in a table
    T(i,1)=M;
end
T=array2table(T,'VariableNames',{'M','p','phi0','GCI12','GCI23','coeff','Check'})
if CFL==0.8
    save('Case1_CFL08');
elseif CFL==0.5
    save('Case2_CFL05');
end
