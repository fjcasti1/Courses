%% Homework 4, Problem 3 - Francisco Castillo'
clear all; close all; clc;
labelfontsize = 14;
linewidth = 2;

N = 200;
dt = 8e-1/N^4;
[D,x] = cheb(N); % D:(N+1)x(N+1), x:(N+1)x1
x = x(2:N);
D2=D^2;
D2 = D2(2:N,2:N);
u = zeros(size(x));
u35 = u;
u0 = u(N/2);
tobs = [3.5 100];
t = 0;
k=1;
igraph = 1000;
tstep = 0;

while u0<5
    if (t+dt>tobs(k))
        dt=tobs(k)-t;
        t=t+dt;
        k=k+1;
    else
        t=t+dt;
        dt = 8e-1/N^4;
    end

    u = u +dt*(D2*u+exp(u));
    u0 = u(N/2);
    tstep = tstep+1;
    if t == tobs(1)
        u35 = u;
    end
    if (mod(tstep,igraph)==0 || round(u0,10)>=5)
        h1 = plot([1;x;-1],[0;u;0],'b*');
        hold on
        plot([1;x;-1],[0;u;0],'b')
        if (u35~=0)
            h2 = plot([1;x;-1],[0;u35;0],'r*');
            plot([1;x;-1],[0;u35;0],'r')
        end
        grid on
        axis([-1 1 0 5])
        xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
        ylabel('$u(x,t)$','interpreter','latex','fontsize',labelfontsize)
        hold off
        shg
    end
end
legend([h1 h2],'$u(x,t_5)$', '$u(x,3.5)$','interpreter','latex')
saveas(gcf,'Latex/FIGURES/P3','png')