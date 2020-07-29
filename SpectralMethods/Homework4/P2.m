%% Homework 4, Problem 2 - Francisco Castillo'
clear all; close all; clc;
labelfontsize = 14;
markersize = 4;
linewidth = 2;

N = 50;
V = [7,8];
levels = log10([1e-2 1e-3 1e-4 1e-5 1e-6]);
% figure('units','normalized','outerposition',[0 0 1 1])
for j=1:length(V)
    v = V(j);
    dt = v/N^2;
    tf = 1;

    [D,x] = cheb(N); % D:(N+1)x(N+1), x:(N+1)x1
    D = D (2:end,2:end);
    u_n = analytic_sol(x(2:end),0);
    u_n1 = analytic_sol(x(2:end),dt);
    u_n2 = analytic_sol(x(2:end),2*dt);

    t = 2*dt;
    while t<tf
        if (t+dt>tf)
            dt=tf-t;
            t=t+dt;
        else
            t=t+dt;
        end

        u_n3 = u_n2+dt*D*(23*u_n2-16*u_n1+5*u_n)/12;
        u = [analytic_sol(1,t);u_n3];

        u_n = u_n1;
        u_n1 = u_n2;
        u_n2 = u_n3;

        res = analytic_sol(x,t);

        subplot(2,2,j)
        plot(x,u,'b*')
        hold on
        plot(x,u,'b')
        plot(x,res,'r*')
        plot(x,res,'r')
        grid on
        axis([-1 1 0 1])
        xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
        ylabel('$u(x,t)$','interpreter','latex','fontsize',labelfontsize)
        hold off
        shg
    end
    %%
    % Adams-Bashforth stability region and pseudo-spectra of this problem
    dt = v/N^2;
    ee = eig(D);
    dtee = dt*ee;
    xr = 6.5*linspace(min(real(dtee)),0.1,50);
    yr = 2*linspace(min(imag(dtee)),max(imag(dtee)),50);
    [xx,yy] = meshgrid(xr,yr);
    zz = xx+1i*yy;
    ps = 0*zz;
    for k=1:numel(zz)
        ps(k) = min(svd(eye(size(D))*zz(k)-dt*D));
    end

    z = exp(1i*pi*(0:200)/100); r = z-1;
    s = (23-16./z+5./z.^2)/12;


    subplot(2,2,j+2)
    contourf(xx,yy,log10(ps),levels)
    c = colorbar;
    hold on
    plot([-8,8],[0,0],'k')
    plot([0,0],[-8,8],'k')
    plot(real(dtee),imag(dtee),'r*','markersize',markersize)
    plot(r./s,'b','linewidth',linewidth)
    xlabel('$\Re(\lambda)$','interpreter','latex','fontsize',labelfontsize)
    ylabel('$\Im(\lambda)$','interpreter','latex','fontsize',labelfontsize)
    axis([-1 0.1 -1 1])
    grid on
end

saveas(gcf,'Latex/FIGURES/P2','png')

function res = analytic_sol(x,t)
    res = zeros(size(x));
    res(x>=-t) = exp(-60*(x(x>=-t)+t-0.5).^2);
    res(x<-t) = exp(-60*(x(x<-t)+t+0.5).^2);
end
