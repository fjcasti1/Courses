%% HOMEWORK 6- Francisco Castillo

clear variables; close all; format long; clc
path='Latex/FIGURES/';

linewidth=1.2;
markersize=7.5;
legendfontsize=12;
axisfontsize=16;

%% Problem 1
% Part a)
dt = 2e-1;
dydt = @(t,y) sin(t)*(y.^2-(cos(t))^2-1);
err=1;
tol=1e-12;
tspan=[0 50];
y0=1;
i=1;
while err>tol
    N=(tspan(2)-tspan(1))/dt(i);
    [t,y] = rk4(dydt,tspan,y0,N);
    err(i)=norm(y-cos(t),inf);
    dt(i+1)=dt(i)/2;
    i=i+1;
end
dt(end)=[];

% Linear fit
P=polyfit(log(dt),log(err),1);
dtfit=2*dt(1):-0.001:dt(end)/2;
errfit=exp(P(2))*dtfit.^(P(1));

figure(1)
loglog(dt,err,'r*','markersize',markersize)
grid on
hold on
loglog(dtfit,errfit,'b-.','linewidth',linewidth)
grid on
xlabel('$\log(\Delta t)$','Interpreter','latex')
ylabel('$\log(\epsilon)$','Interpreter','latex')
set(gca,'fontsize',14)
txt=[path,'problem1_a'];
saveas(gcf,txt,'epsc')

% Part b)
N=2.^(4:12);
for i=1:length(N);
    [t1,y1] = euler(dydt,tspan,y0,N(i));
    [t2,y2] = modEuler(dydt,tspan,y0,N(i));
	[t4,y4] = rk4(dydt,tspan,y0,N(i));
    err1(i)=norm(y1-cos(t1),inf);
    err2(i)=norm(y2-cos(t2),inf);
    err4(i)=norm(y4-cos(t4),inf);
end
figure
loglog(N,err1,2*N,err2,4*N,err4)
grid on
ylabel('$\log(\epsilon)$','Interpreter','latex')
set(gca,'fontsize',12)
xlabel('Function Evaluations','Interpreter','latex')
legend({'Euler','Modified Euler', 'RK4'},'fontsize',14)
txt=[path,'problem1_b'];
saveas(gcf,txt,'epsc')

%% Problem 2
Nx=100; Ny=300;
hx=1/(Nx-1); hy=1/(Ny-1);

D2x=gallery('tridiag',Nx,1,-2,1)/(hx^2);
D2y=gallery('tridiag',Ny,1,-2,1)/(hy^2);

[Vx,LambdaX]=eig(full(D2x));
[Vy,LambdaY]=eig(full(D2y));
[x,y]=meshgrid(linspace(0,1,100),linspace(0,1,300));
f=sin(x).*cos(100*y);
F=Vy'*f*Vx;
U=zeros(size(f));
lambdax=diag(LambdaX);
lambday=diag(LambdaY);
% whos U F lambdax lambday
for j=1:Nx
    U(:,j)=F(:,j)./(lambday+lambdax(j));
end

u=Vy*U*Vx';
% Calculate error
err=norm(u*D2x+D2y*u-f,inf)
figure
surf(x,y,u)
shading interp
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
txt=[path,'problem2_u'];
saveas(gcf,txt,'png')

%% Problem 3
clear variables; close all; clc
DATA=load('velocitydataPlatte.mat');
u=DATA.u;
v=DATA.v;
xx=DATA.xx;
yy=DATA.yy;
fluidflowFDTemperature(xx,yy,u,v)

%% Problem 4
Ra=2e5;
Pr=0.71;
Problem4(Ra,Pr)







