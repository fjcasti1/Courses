%% Homework 4, Problem 4 -  Francisco Castillo
clear all; close all; clc
labelfontsize = 14;
%% 2D wave equation Chebyshev+Leap-frog, ZERO Neumann BC'S
N = 64;
[D,x] = cheb(N);
% x = x(2:end-1);
% D2 = D^2;
% D2 = D2(2:end-1,2:end-1);
 
% 2D grid
[X,Y] = meshgrid(x);
 
% Initial condition
u0 = exp(-40*((X-0.2).^2+Y.^2));
u = u0;
 
h = 1-x(2);
dt = h/2;
t = 0;
tf = 10;
count=0;

while t<tf
    if t+dt>tf
        dt = tf-t;
    else
        dt = h/2;
    end
    uyy = D(:,2:N)*D(2:N,:)*u;
    uxx = u*D(2:N,:)'*D(:,2:N)';
    u2 = 2*u- u0 + dt^2*(uxx+uyy);
    u0 = u;
    u = u2;
    
    if count == 10 
        surf(X,Y,u)
        zlim([-1 1])
        drawnow
        shg
        count = 0;
    end
        
    count = count+1;
    t = t+dt;
    
    if t == tf
        surf(X,Y,u)
        zlim([-1 1])
        drawnow
        shg
        xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
        ylabel('$y$','interpreter','latex','fontsize',labelfontsize)
        zlabel('$u(x,y)$','interpreter','latex','fontsize',labelfontsize)
        set(get(gca,'ZLabel'),'Rotation',0)
        saveas(gcf,'Latex/FIGURES/P4','png')
    end
end
