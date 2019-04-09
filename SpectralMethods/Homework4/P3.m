%% Homework 4, Problem 3 - Francisco Castillo'
clear all; close all; clc;
labelfontsize = 14;
linewidth = 2;

N = 50;

[D,x] = cheb(N); % D:(N+1)x(N+1), x:(N+1)x1
D2=D^2;
D2 = D2(2:N,2:N);
u = zeros(size(x));
u = u(2:N);
uold = u;
tobs = [3.5 5];
t = 0;
k=1;
while t<tobs(end)
    dt = 1/N^3;
    if (t+dt>tobs(k))
        dt=tobs(k)-t;
        t=t+dt;
        k=k+1;
    else
        t=t+dt;
    end
    cccc
    u = uold +dt*(D2*u+exp(u));
    
%     u_n3 = u_n2+dt*D*(23*u_n2-16*u_n1+5*u_n)/12;
%     u = [analytic_sol(1,t);u_n3];
%     u_n = u_n1;
%     u_n1 = u_n2;
%     u_n2 = u_n3;

    uold = u;
    
    plot(x,[0;u;0],'b*')
    hold on
    plot(x,[0;u;0],'b')
    grid on
%     axis([-1 1 0 1])
    xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
    ylabel('$u(x,t)$','interpreter','latex','fontsize',labelfontsize)
    hold off
    shg
    
    cccc
end
xcc
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
contourf(xx,yy,log10(ps))
c = colorbar;
hold on
plot([-8,8],[0,0],'k')
plot([0,0],[-8,8],'k')
plot(real(dtee),imag(dtee),'r*')
plot(r./s,'b','linewidth',linewidth)
xlabel('$\Re(\lambda)$','interpreter','latex','fontsize',labelfontsize)
ylabel('$\Im(\lambda)$','interpreter','latex','fontsize',labelfontsize)
axis([-1 0.1 -1 1])
grid on


saveas(gcf,'Latex/FIGURES/P2','png')