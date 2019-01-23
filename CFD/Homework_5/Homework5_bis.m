%% HOMEWORK 5 - FRANCISCO CASTILLO
%%
%% Defined functions
%
%% Problem 1
format long
clear all; close all; clc
labelfontsize=14;

M=256;
L=2;
alpha = 0.1;
h=L/M;
dt=4*h^2/(2*alpha);
x=linspace(-1,1,M+1)';
g = @(t) 2-sin(3*pi*t/2);
q = @(x,t) 3*pi*sin(pi*sqrt(t)/2).*cos(pi*sqrt(t)/2).*(x.^3-x.^2-x+1)./(2*sqrt(t))...
    +3*pi*cos(3*pi*t/2).*sin(pi*x/2)/2+...
    alpha*(sin(pi*sqrt(t)/2).^2.*(6-18*x)+pi^2*sin(3*pi*t/2).*sin(pi*x/2)/4);

T=2*ones(length(x),1);
Tint=T(2:end-1); %Initialize with initial condition.
xint=x(2:end-1);
qnew=zeros(length(x),1);
time=0;
outputTime=[1 2 3];
endtime=outputTime(end);
j=1;

figure(1)
plot(x,T)
hold on
axis([-1 1 0.5 6.5])
grid on
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
ylabel('$T(x,t)$','fontsize',labelfontsize,'interpreter','latex')
title('Temperature profiles')

while time < endtime
    if (time < outputTime(j) && time+dt >= outputTime(j))
        dt=outputTime(j)-time;
    else
        dt=4*h^2/(2*alpha);
    end
    qold=qnew;
    qnew=q(x,time+dt);
    
    B=alpha*dt/(2*h^2);
    a=-B*ones(length(x)-2,1);
%     a(1)=0;
    b=(1+2*B)*ones(length(x)-2,1);
    b(end)=1+B;
    c=-B*ones(length(x)-2,1);
%     c(end)=0;
    d=zeros(length(x)-2,1);

    for i=2:length(Tint)-1 %interior of the interior for d
        d(i)=Tint(i)+B*(Tint(i+1)-2*Tint(i)+Tint(i-1))+...
           dt*0.5*(qold(i+1)+qnew(i+1));
    end
    d(1)=Tint(1)+B*(Tint(2)-2*Tint(1)+g(time)+g(time+dt))+...
        +dt*0.5*(qold(2)+qnew(2));
    d(end)=Tint(end)+B*(-Tint(end)+Tint(end-1))+...
        +dt*0.5*(qold(end-1)+qnew(end-1));
    %
    Tint(:)=GaussTriSol(a,b,c,d);
    T(1)=g(time+dt);
    T(2:end-1)=Tint(:);
    T(end)=T(end-1);
    if time+dt == outputTime(j)
        j=j+1;
        figure(1)
        plot(x,T)
        axis([-1 1 0.5 6.5])
    end
    time=time+dt;
end
lh = legend({'$T(x,0)$','$T(x,1)$','$T(x,2)$','$T(x,3)$'},...
    'Interpreter','latex','Fontsize',14);
ccc
%% Problem 2
format long
clear all; clc
labelfontsize=14;

M=256;
L=2;
alpha = 0.1;
h=L/M;
dt=4*h^2/(2*alpha);
x=[-1-h/2:h:1+h/2]';
g = @(t) 2-sin(3*pi*t/2);
q = @(x,t) 3*pi*sin(pi*sqrt(t)/2).*cos(pi*sqrt(t)/2).*(x.^3-x.^2-x+1)./(2*sqrt(t))...
    +3*pi*cos(3*pi*t/2).*sin(pi*x/2)/2+...
    alpha*(sin(pi*sqrt(t)/2).^2.*(6-18*x)+pi^2*sin(3*pi*t/2).*sin(pi*x/2)/4);

T=2*ones(length(x),1);
Tint=T(2:end-1); %Initialize with initial condition.
xint=x(2:end-1);
qnew=zeros(length(x),1);
time=0;
outputTime=[1 2 3];
endtime=outputTime(end);
j=1;

figure(2)
plot(x,T)
hold on
axis([-1 1 0.5 6.5])
grid on
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
ylabel('$T(x,t)$','fontsize',labelfontsize,'interpreter','latex')
title('Temperature profiles')

while time < endtime
    if (time < outputTime(j) && time+dt >= outputTime(j))
        dt=outputTime(j)-time;
    else
        dt=4*h^2/(2*alpha);
    end
    qold=qnew;
    qnew=q(x,time+dt);
    
    B=alpha*dt/(2*h^2);
    a=-B*ones(length(x)-2,1);
    a(1)=0;
    b=(1+2*B)*ones(length(x)-2,1);
    b(1)=1+3*B;
    b(end)=1+B;
    c=-B*ones(length(x)-2,1);
    c(end)=0;
    d=zeros(length(x)-2,1);

    for i=2:length(Tint)-1 %interior of the interior for d
        d(i)=Tint(i)+B*(Tint(i+1)-2*Tint(i)+Tint(i-1))+...
           dt*0.5*(qold(i+1)+qnew(i+1));
    end
    d(1)=Tint(1)+B*(Tint(2)-3*Tint(1)+2*g(time)+2*g(time+dt))+...
        +dt*0.5*(qold(2)+qnew(2));
    d(end)=Tint(end)+B*(-Tint(end)+Tint(end-1))+...
        +dt*0.5*(qold(end-1)+qnew(end-1));

    Tint(:)=GaussTriSol(a,b,c,d);
    T(2:end-1)=Tint(:);
    T(1)=2*g(time+dt)-T(2);
    T(end)=T(end-1);
    if time+dt == outputTime(j)
        j=j+1;
        figure(2)
        plot(x,T)
        axis([-1 1 0.5 6.5])
    end
    time=time+dt;
end
lh = legend({'$T(x,0)$','$T(x,1)$','$T(x,2)$','$T(x,3)$'},...
    'Interpreter','latex','Fontsize',14);
% lh = legend({'$T_n(x,0)$','$T_n(x,1)$','$T_n(x,2)$','$T_n(x,3)$',...
%     '$T_c(x,0)$','$T_c(x,1)$','$T_c(x,2)$','$T_c(x,3)$'},...
%     'Interpreter','latex','Fontsize',14);