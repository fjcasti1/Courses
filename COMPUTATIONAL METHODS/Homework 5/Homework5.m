%% HOMEWORK 5 - FCO CASTILLO
clear variables
close all
clc
%% Problem 1a
legendfontsize=14;
axisfontsize=16;

F = @(t,V) [-4*V(2)+V(1).*(1-V(1).^2-V(2).^2);...
                4*V(1)+V(2).*(1-V(1).^2-V(2).^2)];
theta=pi/4;
R=[0.1; 2.4];
tspan=[0 10];

N=300;
a=2.5;
for i=1:length(R)
    V0 = [R(i)*cos(theta); R(i)*sin(theta)];
    [t,V] = rk4(F,tspan,V0,N);
    figure(1)
    l(1)=plot(V(:,1),V(:,2),'r');
    hold on
    l(2)=plot(V0(1),V0(2),'r*');
    axis([-a a -a a])
    grid on
end
% Plot orbits
figure(1)
hold on
l(3)=plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'--b');
xlabel('$x(t)$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
pbaspect([1 1 1])
set(gca,'fontsize',14)
legend([l(1) l(2) l(3)],{' Orbits',' Initial condtions',' Limit circle'},...
 'Interpreter','latex','fontsize',legendfontsize,'location','SouthEast')
txt='Latex/FIGURES/P1_1';
saveas(gcf,txt,'epsc')

% Plot vector field
step=0.2;
[x,y] = meshgrid(-a:step:a,-a:step:a);
u=-4*y+x.*(1-x.^2-y.^2);
v=4*x+y.*(1-x.^2-y.^2);
figure(2)
q=quiver(x,y,u,v,'r');
hold on
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'--b');
xlabel('$x(t)$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
pbaspect([1 1 1])
axis([-a a -a a])
set(gca,'fontsize',14)
set(q,'AutoScale','on', 'AutoScaleFactor', 1.5)
txt='Latex/FIGURES/P1_1field';
saveas(gcf,txt,'epsc')

%% Problem 1b
close all
F = @(t,V) [-4*V(2)+V(1).*(1-V(1).^2-V(2).^2).*(4-V(1).^2-V(2).^2);...
                4*V(1)+V(2).*(1-V(1).^2-V(2).^2).*(4-V(1).^2-V(2).^2)];
theta=pi/5;
R=[0.1; 1.2; 1.8];
tspan=[0 10];
N=300;
a=3;
for i=1:length(R)
    V0 = [R(i)*cos(theta); R(i)*sin(theta)];
    [t,V] = rk4(F,tspan,V0,N);
    figure(1)
    plot(V(:,1),V(:,2),'r')
    hold on
    plot(V0(1),V0(2),'r*')
    axis([-a a -a a])
    grid on
end
Theta=0:pi/4:2*pi;
for i=1:length(Theta)
    V0 = [2*cos(Theta(i)); 2*sin(Theta(i))];
    [t,V] = rk4(F,[0 1],V0,N);
    figure(1)
    l(1)=plot(V(:,1),V(:,2),'r');
    hold on
    l(2)=plot(V0(1),V0(2),'r*');
    axis([-a a -a a])
    grid on
end

figure(1)
hold on
l(3)=plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'--b');
plot(2*cos(0:0.001:2*pi),2*sin(0:0.001:2*pi),'--b');
xlabel('$x(t)$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
pbaspect([1 1 1])
set(gca,'fontsize',14)
%legend([l(1) l(2) l(3)],{' Orbits',' Initial condtions',' Limit circle'},...
% 'Interpreter','latex','fontsize',legendfontsize,'location','SouthEast')
txt='Latex/FIGURES/P1_2';
saveas(gcf,txt,'epsc')

% Plot vector field
step=0.1;
b=2;
[x,y] = meshgrid(-b:step:b,-b:step:b);
u=-4*y+x.*(1-x.^2-y.^2).*(4-x.^2-y.^2);
v=4*x+y.*(1-x.^2-y.^2).*(4-x.^2-y.^2);
figure(2)
q=quiver(x,y,u,v,'r');
set(q,'AutoScale','on', 'AutoScaleFactor', 3)
hold on
% step=0.1;
% b=2.1;
% [x,y] = meshgrid(-b:step:b,-b:step:b);
% u=-4*y+x.*(1-x.^2-y.^2).*(4-x.^2-y.^2);
% v=4*x+y.*(1-x.^2-y.^2).*(4-x.^2-y.^2);
% q=quiver(x,y,u,v,'r');
% step=0.05;
% b=1.2;
% [x,y] = meshgrid(-b:step:b,-b:step:b);
% u=-4*y+x.*(1-x.^2-y.^2).*(4-x.^2-y.^2);
% v=4*x+y.*(1-x.^2-y.^2).*(4-x.^2-y.^2);
% q=quiver(x,y,u,v,'r');
% set(q,'AutoScale','on', 'AutoScaleFactor', 0.5)
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'--b');
plot(2*cos(0:0.001:2*pi),2*sin(0:0.001:2*pi),'--b');
xlabel('$x(t)$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
pbaspect([1 1 1])
grid on
a=2.1;
axis([-a a -a a])
set(gca,'fontsize',14)
txt='Latex/FIGURES/P1_2field';
saveas(gcf,txt,'epsc')

%% Problem 2b
clear variables
close all

N=100;
h=pi/N;

Du = gallery('tridiag',N,-1,1,0)/h; % In sparse form.
Du(:,end)=[];
Dp = gallery('tridiag',N,0,-1,1)/h; % In sparse form.
Dp(end,:)=[];
Z1 = zeros(N-1,N-1);
Z2 = zeros(N,N);

M = [Z1 Dp ; Du Z2];

E=eig(full(M));
E=sort(E);
k=-(N-1):(N-1);
E_analytic=2*1i/h*sin(k*h/2)';
E_analytic=sort(E_analytic); 

difference=norm(E-E_analytic,inf) 
% The expression in part b is correct
%% Problem 2e and 2f
a=-20;
b=20;
c=-20;
d=20;
dx=0.1;
dy=dx;

[zr,zi]=meshgrid(a:dx:b,c:dy:d);
z=zr+1i*zi;
g=1+1/6*(6*z+3*z.^2+z.^3+z.^4/4);
dt=1.405*h;
lambda=dt*E;

figure
contourf(zr,zi,abs(g),[-inf 1 inf], 'k'), colorbar
hold on
grid on
plot(real(lambda),imag(lambda),'r*')
xlabel('$x(t)$','interpreter','latex','fontsize',16)
ylabel('$y(t)$','interpreter','latex','fontsize',16)
axis('image', [-4 1 -3 3])
set(gca,'fontsize',12)
txt='Latex/FIGURES/P2_ef';
saveas(gcf,txt,'epsc')
%% Problem 2g
clear variables
close all
format long
clc

linewidth=1.7;
legendfontsize=14;
axisfontsize=14;

N=100;
L=pi;
h=L/N;
xu=linspace(0,pi,N+1);
xp=linspace(0+h/2,pi-h/2,N);

Du = gallery('tridiag',N,-1,1,0); % In sparse form.
Du(:,end)=[];
Du=Du/h;
Dp = gallery('tridiag',N,0,-1,1); % In sparse form.
Dp(end,:)=[];
Dp=Dp/h;
Z1 = zeros(N-1,N-1);
Z2 = zeros(N,N);
M = [Z1 Dp ; Du Z2];

u0 = exp(-30*(xu-pi/2).^2)';
p0 = zeros(N,1);
V0 = [u0(2:end-1) ; p0];

t = 0:pi/200:2*pi;
dydt = @(t,V) M*V;
[t,V] = rk4(dydt,[0 2*pi],V0,length(t)-1);
% Plots
figure(1)
plot(xu,[0 V(101,1:N-1) 0],'linewidth',linewidth)
hold on
plot(xp,V(101,N:2*N-1),'linewidth',linewidth)
grid on
axis([0 pi -1 1])
xlabel('$x(t)$','interpreter','latex')
legend({' $u(x,t)$',' $p(x,t)$'},...
 'Interpreter','latex','fontsize',legendfontsize)
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P2_g1';
saveas(gcf,txt,'epsc')


figure(2)
plot(xu,[0 V(201,1:N-1) 0],'linewidth',linewidth)
hold on
plot(xp,V(201,N:2*N-1),'linewidth',linewidth)
grid on
axis([0 pi -1 1])
xlabel('$x(t)$','interpreter','latex')
legend({' $u(x,t)$',' $p(x,t)$'},...
 'Interpreter','latex','fontsize',legendfontsize)
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P2_g2';
saveas(gcf,txt,'epsc')

figure(3)
plot(xu,[0 V(301,1:N-1) 0],'linewidth',linewidth)
hold on
plot(xp,V(301,N:2*N-1),'linewidth',linewidth)
grid on
axis([0 pi -1 1])
xlabel('$x(t)$','interpreter','latex')
legend({' $u(x,t)$',' $p(x,t)$'},...
 'Interpreter','latex','fontsize',legendfontsize)
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P2_g3';
saveas(gcf,txt,'epsc')

figure(4)
plot(xu,[0 V(401,1:N-1) 0],'linewidth',linewidth)
hold on
plot(xp,V(401,N:2*N-1),'linewidth',linewidth)
grid on
axis([0 pi -1 1])
xlabel('$x(t)$','interpreter','latex')
legend({' $u(x,t)$',' $p(x,t)$'},...
 'Interpreter','latex','fontsize',legendfontsize)
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P2_g4';
saveas(gcf,txt,'epsc')
%% Problem 2h
t = 0:1.48*h:2*pi;
dydt = @(t,V) M*V;
[t,V] = rk4(dydt,[0 2*pi],V0,length(t)-1);
figure 
plot(xu,[0 V(end,1:N-1) 0])
hold on
plot(xp,V(end,N:2*N-1))
grid on
% axis([0 pi -1 1])
xlim([0 pi])
xlabel('$x(t)$','interpreter','latex')
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P2_h';
saveas(gcf,txt,'epsc')
%% Problem 2i
clear variables
close all
format long
clc
%Calculate exact solution with a large number of N
axisfontsize=14;
L=pi;

% Check order of convergence
Nvector=[50 100 200 400 800 1600];
for i=1:length(Nvector)
    N=Nvector(i)
    h=L/N;
    xu=linspace(0,pi,N+1);
    xp=linspace(0+h/2,pi-h/2,N);
    
    Du = gallery('tridiag',N,-1,1,0); % In sparse form.
    Du(:,end)=[];
    Du=Du/h;
    Dp = gallery('tridiag',N,0,-1,1); % In sparse form.
    Dp(end,:)=[];
    Dp=Dp/h;
    Z1 = zeros(N-1,N-1);
    Z2 = zeros(N,N);
    M = [Z1 Dp ; Du Z2];
    
    u0 = exp(-30*(xu-pi/2).^2)';
    p0 = zeros(N,1);
    V0 = [u0(2:end-1) ; p0];
    dt=h;
    t = 0:dt:2*pi;
    dydt = @(t,V) M*V;
    [t,V] = rk4(dydt,[0 2*pi],V0,length(t)-1);
    
    difNormu(i)=norm(V(end,1:N-1)'-u0(2:end-1),inf);
    index=find(xp>=pi/2,1);
    difNormp(i)=norm(V(end,N:2*N-1)'-p0,inf);
end
% Order of convergence of u and p
hvector=pi/1600:1e-5:pi/50;
NvectorCont=20:500:7e3;

% Fitting curves for the error.
% Ignore first 2 points. No asymptotic regime.

cu=polyfit(log(Nvector(3:end)),log(difNormu(3:end)),1);
ufit=exp(cu(2))*NvectorCont.^cu(1);
cp=polyfit(log(Nvector(3:end)),log(difNormp(3:end)),1);
pfit=exp(cp(2))*NvectorCont.^cp(1);
cu(1)
cp(1)

figure
loglog(Nvector,difNormu,'r*')
hold on
loglog(NvectorCont,ufit,'linewidth',2)
grid on
xlim([30 2000])
set(gca,'fontsize',axisfontsize)
xlabel('$N$','interpreter','latex','fontsize',16)
ylabel('$e_u$','interpreter','latex','fontsize',16)
txt='Latex/FIGURES/P2_iu';
saveas(gcf,txt,'epsc')

figure
loglog(Nvector,difNormp,'r*')
hold on
loglog(NvectorCont,pfit,'linewidth',2)
grid on
xlim([30 2000])
set(gca,'fontsize',axisfontsize)
xlabel('$N$','interpreter','latex','fontsize',16)
ylabel('$e_p$','interpreter','latex','fontsize',16)
txt='Latex/FIGURES/P2_ip';
saveas(gcf,txt,'epsc')
%% Problem 2j
clear variables
close all
clc

legendfontsize=14;
axisfontsize=14;

L=pi;
N=100;
h=L/N;

theta=linspace(0,2*pi,N);
g=exp(1i*theta);
z=12*(g.^3-g.^2)./(23*g.^2-16*g+5);

xu=linspace(0,pi,N+1);
xp=linspace(0+h/2,pi-h/2,N);

Du = gallery('tridiag',N,-1,1,0); % In sparse form.
Du(:,end)=[];
Du=Du/h;
Dp = gallery('tridiag',N,0,-1,1); % In sparse form.
Dp(end,:)=[];
Dp=Dp/h;
Z1 = zeros(N-1,N-1);
Z2 = zeros(N,N);
M = [Z1 Dp ; Du Z2];

e=eig(full(M));
dt=0.36*h;
lambdadt=dt*e;

%Plot eigenvalues in stability region
figure
plot(z,'b','linewidth',2)
hold on
plot(real(lambdadt),imag(lambdadt),'r*')
grid on
xlabel('$\Re(z)$','interpreter','latex'...
    ,'fontsize',16)
ylabel('$\Im(z)$','interpreter','latex'...
    ,'fontsize',16)
txt='Latex/FIGURES/P2_jRegion';
saveas(gcf,txt,'epsc')


u0 = exp(-30*(xu-pi/2).^2)';
p0 = zeros(N,1);
V0 = [u0(2:end-1) ; p0];


F2=@(t,v) M*v;

[t,V]=rk4(F2,[0 2*dt],V0,2);

F1=@(v) M*v;
f0=F1(V(1,:)');
f1=F1(V(2,:)');
f=F1(V(3,:)');

V=V(3,:)'+dt/12*(23*f-16*f1+5*f0);

t=t(end);
outputTime=[pi/2 pi 3*pi/2 2*pi];
n=1;
while t<2*pi
    if (t < outputTime(n) && t+dt >= outputTime(n))
        dt=outputTime(n)-t;
    else
        dt=0.36*h;
    end
    f=F1(V);
    V=V+dt/12*(23*f-16*f1+5*f0);
    %Advance variables
    f0=f1;
    f1=f;
    t=t+dt;
    if (t==outputTime(n))
        figure
        plot(xu',[0; V(1:N-1); 0],'linewidth',2)
        hold on
        plot(xp',V(N:2*N-1),'linewidth',2)
        grid on
        legend({' $u(x,t)$',' $p(x,t)$'},...
        'Interpreter','latex','fontsize',legendfontsize)
        set(gca,'fontsize',axisfontsize)
        xlabel('$x(t)$','interpreter','latex','fontsize',16)
        axis([0 pi -1 1])
        txt=['Latex/FIGURES/P2_j' num2str(n)];
        saveas(gcf,txt,'epsc')
        n=n+1;
    end
end

%% Problem 3
clear variables
close all
clc
N=300;
theta = linspace(0,2*pi,N);
g=exp(1i*theta);
a=115*g.^2-80*g+25;
b=13*g.^2-g;
c=g.^2-g.^3;
sq=sqrt(b.^2-4*a.*c);
zp=12*(-b+sq)./(2*a);
zn=12*(-b-sq)./(2*a);
figure
plot(real(zp),imag(zp),'b.')
hold on
plot(real(zn),imag(zn),'b.')
grid on
axis('image', [-2.5 0.5 -1.5 1.5])
xlabel('$\Re(z)$','interpreter','latex'...
    ,'fontsize',16)
ylabel('$\Im(z)$','interpreter','latex'...
    ,'fontsize',16)
txt='Latex/FIGURES/P3';
saveas(gcf,txt,'epsc')