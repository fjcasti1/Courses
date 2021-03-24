clear all; close all; clc

%% Parameter values
v = -1;
d = 0.01;
L = 1;

T = 0.5;
Ntsteps = 100;
dt = T/Ntsteps;
t = linspace(0,T,Ntsteps+1);

R = 10000; % Number of realizations

M = 100; % Number of grid cells
dx = 1/M;
x = linspace(0,1,(M+1))';

% Distribution size
b = sqrt(24*d/dt);
a = v - 0.5*b;

%% Initial condition
xi_0 = 0.5;
u_0  = calculate_density(x,xi_0);
u(:,1) = u_0;

for r = 1:R
    
    n = 1; % time-step counter
    xi(1,r) = xi_0;
    u(:,1,r) = u_0;
    
    while t(n)<T
        eta = rand(1);
        q = a + eta*b;
        xi_new = xi(n,r) + dt*q;
        % Periodic Boundary conditions
        if (xi_new > L)
            xi(n+1,r) = xi_new - L;
        elseif (xi_new < 0)
            xi(n+1,r) = xi_new + L;
        else
            xi(n+1,r) = xi_new;
        end
        u(:,n+1,r) = calculate_density(x,xi(n+1,r));
        
        n = n + 1;
    end
end

% Average across realizations
u = mean(u,3);

% Plot solution
[tt,xx] = meshgrid(t,x);
figure(2)
plot3(tt,xx,u)
saveas(gcf,'../figures/p2_sol','png')

function u = calculate_density(x, xi)
    M = length(x)-1;
    dx = (x(end)-x(1))/M;
    x_staggered = linspace(x(1)-dx,x(end)+dx,M+2);
    N = 1; % Only one particle moving
%     dx = 1;
    for m=1:M+1
        I=find((x_staggered(m)<xi)&(xi<x_staggered(m+1)));
        u(m)=length(I)/N/dx;
    end
end