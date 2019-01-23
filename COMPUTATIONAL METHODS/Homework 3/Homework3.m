%% Problem 1
clear all; close all; format long; clc
legendfontsize=14;
axisfontsize=16;
f = @(x) (exp(3*x).*(sin(200*x.^2)))./(1 + 20*x.^2);
xx=0:0.001:1;
for j=4:14
    n=2^j;
    x=linspace(0,1,n+1);
    fx=f(x);
    fxx=spline(x,fx,xx);
    error=max(abs(f(xx)-fxx));
    
    figure(1)
    loglog(n,error,'r*')
    hold on
    grid on
    if j==14
        xlabel('$n$','Interpreter','latex')
        set(gca,'fontsize',14)
        legend({'Error'},...
         'Interpreter','latex','fontsize',legendfontsize)
        txt='Latex/FIGURES/P1_Error';
        saveas(gcf,txt,'epsc')
    end
    
    figure(2)
    semilogy(n,error,'r*')
    hold on
    grid on
    if j==14
        xlabel('$n$','Interpreter','latex')
        set(gca,'fontsize',14)
        legend({'Error'},...
         'Interpreter','latex','fontsize',legendfontsize)
        txt='Latex/FIGURES/P1_Error_semilog';
        saveas(gcf,txt,'epsc')
    end
    
    figure
    plot(xx,f(xx))
    hold on
    plot(xx,fxx)
    grid on
    xlabel('$x$','Interpreter','latex')
    legend({'$f(x)$','Spline'},...
        'Interpreter','latex','fontsize',legendfontsize)
    set(gca,'fontsize',axisfontsize)
    txt=['Latex/FIGURES/P1_n=',num2str(n)];
    saveas(gcf,txt,'epsc')
end

%% Problem 2
clear all; close all
legendfontsize=14;
axisfontsize=14;
labelfontsize=16;
%% Original code, non-sparse matrices.
time=zeros(8,1);
for j=1:8
    tic
    N = 100*j;
    x = linspace(-1,1,N)';
    dx = x(2)-x(1);

    unos = ones(N-2,1);
    D2 = ( diag(unos(1:end-1),-1)+diag(unos(1:end-1),1)-2*diag(unos,0) )/dx^2;
    B = zeros(N-2,2);
    B(1,1) = 1/dx^2; B(end,end) =1/dx^2;
    % initial condition
    u0 = sin(pi*x)+0.2*sin(5*pi*x);
    % boundary conditions
    g1 = @(t) .5*sin(50*t);
    g2 = @(t) 0*t;

    t = 0:0.001:.5;

    [T,U] = ode45(@(t,u) D2*u+B*[g2(t);g1(t)], t, u0(2:end-1));
%     for k = 2:length(t)
%         plot(x,[0;U(k,:)';g1(t(k))],'*-')
%         ylim([-1.5 1.5])
%         shg
%         drawnow
%     end
    time(j)=toc;
    if N==100
        figure
        [T,X]= meshgrid(t,x(2:end-1));
        surf(T,X,U','edgecolor','none')
        xlabel('$x$','fontsize',labelfontsize,...
            'interpreter','latex')
        ylabel('$y$','fontsize',labelfontsize,...
            'interpreter','latex')
        zlabel('$\phi_V(x,y)$','fontsize',labelfontsize,...
            'interpreter','latex')
        set(gca,'fontsize',axisfontsize)
        txt='Latex/FIGURES/P2_surf1';
        saveas(gcf,txt,'epsc')
    end
end

%% Modified code, sparse matrices.
sparsetime=zeros(8,1);
for j=1:8
    tic
    N = 100*j;
    x = linspace(-1,1,N)';
    dx = x(2)-x(1);

    unos = ones(N-2,1);
    D2 = sparse(( diag(unos(1:end-1),-1)+diag(unos(1:end-1),1)-2*diag(unos,0) )/dx^2);
    B = sparse(zeros(N-2,2));
    B(1,1) = 1/dx^2; B(end,end) =1/dx^2;
    % initial condition
    u0 = sin(pi*x)+0.2*sin(5*pi*x);
    % boundary conditions
    g1 = @(t) .5*sin(50*t);
    g2 = @(t) 0*t;

    t = 0:0.001:.5;

    [T,U] = ode45(@(t,u) D2*u+B*[g2(t);g1(t)], t, u0(2:end-1));
%     for k = 2:length(t)
%         plot(x,[0;U(k,:)';g1(t(k))],'*-')
%         ylim([-1.5 1.5])
%         shg
%         drawnow
%     end
    sparsetime(j)=toc;
    if N==100
        figure
        [T,X]= meshgrid(t,x(2:end-1));
        surf(T,X,U','edgecolor','none')
        txt='Latex/FIGURES/P2_surf2';
        saveas(gcf,txt,'epsc')
    end
end
%%
figure
Nj=100:100:800;
plot(Nj,time,'r*',Nj,sparsetime,'b*')
xlabel('$N$','interpreter','latex')
set(gca,'fontsize',axisfontsize)
legend('Non-Sparse','Sparse')
grid on
txt='Latex/FIGURES/P2_times';
saveas(gcf,txt,'epsc')

%% Problem 3
% Simple code to solve the diffusion equation
% u_t = u_xx -1<x<1, with Neumann boundary conditions
clear all; close all; clc
legendfontsize=14;
axisfontsize=14;
labelfontsize=16;
N = 50;
L=2;
h=L/N;
x=-1:h:1;
D2 = gallery('tridiag',N+1,1,-2,1); % In sparse form.
%%

% initial condition
u0 = cos(pi*x)+0.2*cos(5*pi*x)+1;
% boundary conditions
g1 = @(t) 0*t+1;
g2 = @(t) 0*t+1;

A = D2(2:N,2:N)/h^2;
A(1,1)=-1/h^2;
A(end,end)=-1/h^2;
B = zeros(N-1,2);
B(1,1) = -1/h;
B(end,end) = 1/h;

t = 0:0.001:2;
%%
[T,U] = ode45(@(t,u) A*u+B*[g1(t);g2(t)], t, u0(2:end-1));
% for k = 2:length(t)
%     plot(x(2:end-1),U(k,:)','*-')
%     ylim([-1 2])
%     shg
%     drawnow
% end
[T,X]= meshgrid(t,x(2:end-1));
surf(T,X,U','edgecolor','none')
xlabel('$t$','fontsize',labelfontsize,...
            'interpreter','latex')
ylabel('$x$','fontsize',labelfontsize,...
            'interpreter','latex')
zlabel('$u(x,t)$','fontsize',labelfontsize,...
            'interpreter','latex')
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P3_surf';
saveas(gcf,txt,'epsc')
%% Problem 4
clear all;close all;format long;clc
legendfontsize=14;
axisfontsize=14;
labelfontsize=16;
N = 80;
L=2;
h=L/N;
xu=linspace(-1,1,N+1);
xp=linspace(-1+h/2,1-h/2,N);

Du = gallery('tridiag',N,-1,1,0); % In sparse form.
Du(:,end)=[];
Du=Du/h;
Dp = gallery('tridiag',N,0,-1,1); % In sparse form.
Dp(end,:)=[];
Dp=Dp/h;
Z1 = zeros(N-1,N-1);
Z2 = zeros(N,N);
M = [Z1 Dp ; Du Z2];

% initial condition
u0 = exp(-32*xu.^2)';
p0 = zeros(N,1);
v0 = [u0(2:end-1) ; p0];

t = 0:0.01:4;

[T,V] = ode45(@(t,v) M*v, t, v0);
for k = 2:length(t)
    plot(xu(2:end-1),V(k,1:N-1)','b*-',xp,V(k,N:2*N-1)','r*-')
    grid on
    ylim([-1 1])
    shg
    drawnow
end
%%

[T,X] = meshgrid(t,xu(2:end-1));
figure
contourf(T,X,round(V(:,2:N)',3))
xlabel('$t$','fontsize',labelfontsize,...
            'interpreter','latex')
ylabel('$x$','fontsize',labelfontsize,...
            'interpreter','latex')
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P4_contour1';
saveas(gcf,txt,'epsc')

[T,X] = meshgrid(t,xp);
figure
contourf(T,X,round(V(:,N:end)',3))
xlabel('$t$','fontsize',labelfontsize,...
            'interpreter','latex')
ylabel('$x$','fontsize',labelfontsize,...
            'interpreter','latex')
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P4_contour2';
saveas(gcf,txt,'epsc')
%%
figure
plot(xu(2:end-1),V(find(t==0.5),1:N-1),'b*')
hold on
plot(xp,V(find(t==0.5),N:2*N-1),'r*')
grid on
axis([-1 1 -1 1])
xlabel('$x$','fontsize',labelfontsize,...
            'interpreter','latex')
legend({'$u(0.5,x)$','$p(0.5,x)$'},...
        'Interpreter','latex','fontsize',legendfontsize)
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P4_05';
saveas(gcf,txt,'epsc')

figure
plot(xu(2:end-1),V(find(t==1),1:N-1),'b*')
hold on
plot(xp,V(find(t==1),N:2*N-1),'r*')
grid on
axis([-1 1 -1 1])
xlabel('$x$','fontsize',labelfontsize,...
            'interpreter','latex')
legend({'$u(1,x)$','$p(1,x)$'},...
        'Interpreter','latex','fontsize',legendfontsize)
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P4_1';
saveas(gcf,txt,'epsc')

figure
plot(xu(2:end-1),V(find(t==1.5),1:N-1),'b*')
hold on
plot(xp,V(find(t==1.5),N:2*N-1),'r*')
grid on
axis([-1 1 -1 1])
xlabel('$x$','fontsize',labelfontsize,...
            'interpreter','latex')
legend({'$u(1.5,x)$','$p(1.5,x)$'},...
        'Interpreter','latex','fontsize',legendfontsize)
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P4_15';
saveas(gcf,txt,'epsc')

figure
plot(xu(2:end-1),V(find(t==2),1:N-1),'b*')
hold on
plot(xp,V(find(t==2),N:2*N-1),'r*')
grid on
axis([-1 1 -1 1])
xlabel('$x$','fontsize',labelfontsize,...
            'interpreter','latex')
legend({'$u(2,x)$','$p(2,x)$'},...
        'Interpreter','latex','fontsize',legendfontsize)
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P4_2';
saveas(gcf,txt,'epsc')


u4=V(length(t),1:N-1)';
figure
plot(xu(2:end-1),u0(2:end-1),'b',xu(2:end-1),u4,'r*')
grid on
axis([-1 1 -0.3 1.3])
xlabel('$x$','fontsize',labelfontsize,...
            'interpreter','latex')
legend({'$u(0,x)$','$u(4,x)$'},...
        'Interpreter','latex','fontsize',legendfontsize)
set(gca,'fontsize',axisfontsize)
txt='Latex/FIGURES/P4_periodic';
saveas(gcf,txt,'epsc')
err=norm(u0(2:end-1)-u4,inf)
