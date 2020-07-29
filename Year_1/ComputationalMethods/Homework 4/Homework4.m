%% HOMEWORK 4 - FRANCISCO CASTILLO
clear all;close all;clc
format long
labelfontsize=20;
axisfontsize=12;
%% Problem 1
xx = linspace(-3,3,1e4);
f = @(x) (x - 1).^2.*exp(x);
df = @(x) (x.^2 - 1).*exp(x);

tic
[x,niter,conv]=newton(f,df,2,1e-12,100)
time=toc
% To obtain the roots, modified newton function
[x,niter,conv]=newton_vector(f,df,2,1e-12,100);

figure
plot(xx,f(xx),'linewidth',1.5)
hold on
plot(x,0,'r*')
grid on
axis([-2 2.5 -1 7.5])
set(gca,'fontsize',axisfontsize)
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
ylabel('$y$','interpreter','latex','fontsize',labelfontsize)
saveas(gcf,'Latex/FIGURES/prob1roots','epsc')


%% Problem 2
n = 200;

[x,y]=meshgrid(linspace(-3,3,n));
z = x+y*1i;
z_roots = [1, -1, -1i, 1i, ...
    sqrt(2)*(1+1i),sqrt(2)*(1-1i),sqrt(2)*(-1+1i),sqrt(2)*(-1-1i)];

% We now compute the roots using Newton's method:
r = nan(size(z)); niter = r; conv = r; % allocate memory
f = @(z) z^8+15*z^4-16;
df = @(z) 8*z^7+60*z^3;
for k = 1:n
    for j = 1:n
        [r(k,j),niter(k,j),conv(k,j)] = newton(f,df,z(k,j),1e-16,200);
    end
end
%%
% Plotting the angle of the Newton's method solution, we obtain a fractal:
close all
figure
pcolor(x, y, round(angle(r),3));
mp = colormap;
colormap jet
shading interp
hold on
plot(z_roots, 'k.','markersize',20); axis equal; axis tight; hold off;
axis([-2 2 -2 2])
xlabel('$\Re(x_k)$','interpreter','latex','fontsize',labelfontsize)
ylabel('$\Im(x_k)$','interpreter','latex','fontsize',labelfontsize)
saveas(gcf,'Latex/FIGURES/colormap','epsc')

% Plot the number of iterations required in a colormap
figure
pcolor(x, y, niter);
shading flat, hold on
plot(z_roots, 'k.','markersize',20); axis equal; axis tight; hold off;
axis([-2 2 -2 2])
colorbar
xlabel('$\Re(x_k)$','interpreter','latex','fontsize',labelfontsize)
ylabel('$\Im(x_k)$','interpreter','latex','fontsize',labelfontsize)
saveas(gcf,'Latex/FIGURES/colormap2','epsc')

% Plot the roots in the complex plane
figure
plot(real(z_roots),imag(z_roots),'b.','markerSize',12)
grid on
axis([-1.5 1.5 -1.5 1.5])
set(gca,'fontsize',axisfontsize)
xlabel('$\Re(x_k)$','interpreter','latex','fontsize',labelfontsize)
ylabel('$\Im(x_k)$','interpreter','latex','fontsize',labelfontsize)
saveas(gcf,'Latex/FIGURES/prob2roots','epsc')
%% Problem 3
close all
clc
format long
% Find an answer with sufficient level of convergence
tol=1e-4;
i=1;
N=6;
err=1;
yy=zeros(6+1,1);
while err>tol
    i=i+1;
    y_old=yy(N/6);
    N = i*6;
    xi=linspace(0,pi,N+1)';
    [F,J,y0]=FandJ(N);
    [yy,niter] = newton(F,J,y0,1e-14);
    y_new=yy(N/6);
    % Check solution
    err=max(abs(y_new-y_old));
end
% Results
nodes=N+1
err
% Run it for N=60
N=60;
xi=linspace(0,pi,N+1)';
[F,J,y0]=FandJ(N);
[yy,niter] = newton(F,J,y0,1e-14);
y_new=yy(N/6);
% Plot the solution
figure
plot(xi,[0; yy; 1],'linewidth',1.5)
hold on
plot(xi(N/6+1),y_new,'r.','markerSize',12)
grid on
leg1= legend('Solution', 'Value at $\pi/6$');
set(leg1,'interpreter','latex','FontSize',17);
set(gca,'fontsize',axisfontsize)
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
ylabel('$y$','interpreter','latex','fontsize',labelfontsize)
saveas(gcf,'Latex/FIGURES/prob3sol','epsc')

%% Problem 4
close all, clc
labelfontsize=17;
axisfontsize=14;
sigma = 10;
beta = 8/3;
rho = 28;
phi0=[3 3 20];
eps = 0.001;
tol = odeset('RelTol',1e-6,'AbsTol',1e-6);
f = @(t,a) [-sigma*a(1) + sigma*a(2);...
    rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
%%% ODE 45 %%%
% Plot with tolerence 1e-6, ODE45
[t,phi45] = ode45(f,[0 100],phi0,tol);     % Runge-Kutta 4th/5th order ODE solver
figure
plot3(phi45(:,1),phi45(:,2),phi45(:,3))
grid on
view([50 20])
set(gca,'fontsize',axisfontsize)
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
ylabel('$y$','interpreter','latex','fontsize',labelfontsize)
zlabel('$z$','interpreter','latex','fontsize',labelfontsize)
set(gca,'fontsize',axisfontsize)
set(get(gca,'ZLabel'),'Rotation',10)
saveas(gcf,'Latex/FIGURES/prob4_ode45atractor','epsc')
% Perturb the initial conditions, ODE45
[t,phi45b] = ode45(f,[0 100],[3+eps 3+eps 20+eps],tol);     % Runge-Kutta 4th/5th order ODE solver
errp_ode45 = norm(f(100,phi45) - f(100,phi45b), inf)

%%% ODE 113 %%%
% Plot with tolerence 1e-6, ODE113
[t,phi113] = ode113(f,[0 100],phi0,tol);     % Runge-Kutta 4th/5th order ODE solver
figure
plot3(phi113(:,1),phi113(:,2),phi113(:,3))
grid on
view([50 20])
set(gca,'fontsize',axisfontsize)
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
ylabel('$y$','interpreter','latex','fontsize',labelfontsize)
zlabel('$z$','interpreter','latex','fontsize',labelfontsize)
set(gca,'fontsize',axisfontsize)
set(get(gca,'ZLabel'),'Rotation',10)
saveas(gcf,'Latex/FIGURES/prob4_ode113atractor','epsc')
% Perturb the initial conditions, ODE115
[t,phi113b] = ode113(f,[0 100],[3+eps 3+eps 20+eps],tol);     % Runge-Kutta 4th/5th order ODE solver
errp_ode113 = norm(f(100,phi113) - f(100,phi113b), inf)

%%% ODE 15s %%%
% Plot with tolerence 1e-6, ODE15s
[t,phi15] = ode15s(f,[0 100],phi0,tol);     % Runge-Kutta 4th/5th order ODE solver
figure
plot3(phi15(:,1),phi15(:,2),phi15(:,3))
grid on
view([50 20])
set(gca,'fontsize',axisfontsize)
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
ylabel('$y$','interpreter','latex','fontsize',labelfontsize)
zlabel('$z$','interpreter','latex','fontsize',labelfontsize)
set(gca,'fontsize',axisfontsize)
set(get(gca,'ZLabel'),'Rotation',10)
saveas(gcf,'Latex/FIGURES/prob4_ode15atractor','epsc')
% Perturb the initial conditions, ODE15s
[t,phi15b] = ode15s(f,[0 100],[3+eps 3+eps 20+eps],tol);     % Runge-Kutta 4th/5th order ODE solver
errp_ode15 = norm(f(100,phi15) - f(100,phi15b), inf)
%
% Plot all together
%figure('units','normalized','outerposition',[0 0 1 1])
figure
plot3(phi45(:,1),phi45(:,2),phi45(:,3))
hold on
plot3(phi113(:,1),phi113(:,2),phi113(:,3))
plot3(phi15(:,1),phi15(:,2),phi15(:,3))
grid on
view([50 20])
set(gca,'fontsize',axisfontsize)
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
ylabel('$y$','interpreter','latex','fontsize',labelfontsize)
zlabel('$z$','interpreter','latex','fontsize',labelfontsize)
set(get(gca,'ZLabel'),'Rotation',10)
leg1=legend({'ode45','ode113','ode15s'},'location','northeast');
set(leg1,'interpreter','latex','FontSize',17);
saveas(gcf,'Latex/FIGURES/prob4differences','epsc')

%% Problem 5
n=10;
j=0:n;
x=cos(j*pi/n);
w_a=zeros(n+1,1);
w_b=zeros(n+1,1);
for j=0:n
    xj=x(j+1);
    xk=x;
    xk(j+1)=[];
    w_a(j+1)=1/prod(xj-xk);
end
w_b(1)=2^(n-1)/(2*n);
for j=2:n
    w_b(j)=2^(n-1)*(-1)^(j-1)/n;
end
w_b(n+1)=2^(n-1)*(-1)^n/(2*n);
err=norm(w_a-w_b,inf)