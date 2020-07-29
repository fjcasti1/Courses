%% APM 506 HOMEWORK 1 - FRANCISCO CASTILLO
clear all; close all; clc
format long

%% Problem 1
intf=zeros(1,21);
for k=-10:1:10
    f = chebfun(@(x) exp(-1i*k*x)*exp(sin(x)),[-pi,pi]);
    intf(k+11)=integral(f,-pi,pi);
end
intf'

%% Problem 2
x=-1:0.005:1;
x_0=0.4;
f = chebfun(@(x) exp(x)*sin(3*x));
df = chebfun(@(x) exp(x)*(sin(3*x)+3*cos(3*x)));
ddf = chebfun(@(x) exp(x)*(6*cos(3*x)-8*sin(3*x)));

H=[0.1;0.05;0.025];
ddf_approx=zeros(3,1);
for k=1:3
    h=H(k);
    ddf_approx(k)=(15*(f(x_0-h)+f(x_0+h))-1.5*(f(x_0-2*h)+f(x_0+2*h))+(1/9)*(f(x_0-3*h)+f(x_0+3*h))-(245/9)*f(x_0))/(10*h^2);
end
% figure(1)
% figure('units','normalized','outerposition',[0 0 1 1])
figure(1)
plot(x,f(x),x,df(x),x,ddf(x),x_0,ddf(x_0),'*',x_0,ddf_approx(1),'s',x_0,ddf_approx(2),'s',x_0,ddf_approx(3),'s')
grid on
legend({'$f(x)$','$f''(x)$','$f''''(x)$'},'Interpreter','latex')
xlabel('$x$','Interpreter','latex')
set(gca,'fontsize',18)
saveas(gcf,'IMAGES/problem2_1','epsc')

P=polyfit(log(H),log(abs(ddf(x_0)-ddf_approx)),1);
h=0.01:1e-3:0.15;
figure(2)
loglog(H,abs(ddf(x_0)-ddf_approx),'*')
grid on
xlabel('$h$ (log scale)','fontsize',20,'interpreter','latex')
ylabel('Error (log scale)','fontsize',20,'interpreter','latex')
set(gca,'fontsize',14)
hold on
loglog(h,exp(P(1)*log(h)+P(2)))
grid on
axis([2e-2 1.2e-1 1e-12 1e-2])
saveas(gcf,'IMAGES/problem2_2','epsc')
P
ccc
%% Problem 3
clear all
close all
format long
x_0=pi/3;
f = @(x) tan(x);
fp = @(x) (1/cos(x))^2;
H=logspace(-1,-14,100);
for k = 1:length(H)
    h = H(k);
    df1(k) = (f(x_0+h)-f(x_0))/h;
    df2(k) = (f(x_0+h)-f(x_0-h))/(2*h);
    df4(k) = (-5/60*f(x_0+2*h)+2/3*f(x_0+h)-2/3*f(x_0-h)+5/60*f(x_0-2*h))/h;
    df6(k) = 0.1*(15*(f(x_0+h)-f(x_0-h))/(2*h)-6*(f(x_0+2*h)-f(x_0-2*h))/(4*h)+(f(x_0+3*h)-f(x_0-3*h))/(6*h));
end
figure('units','normalized','outerposition',[0 0 1 1])
loglog(H,abs(df1-fp(x_0)),'*',H,abs(df2-fp(x_0)),'*',H,abs(df4-fp(x_0)),'*',H,abs(df6-fp(x_0)),'*',H,H/2,'--',H,H.^(2)/6,'--',H,H.^(4)/30,'--',H,H.^6/140,'--')
set(gca,'fontsize',14)
grid on
xlabel('$h$ (log scale)','fontsize',20,'interpreter','latex')
ylabel('Error (log scale)','fontsize',20,'interpreter','latex')
saveas(gcf,'IMAGES/problem3','epsc')
h1=H(find(abs(df1-fp(x_0))==min(abs(df1-fp(x_0)))))
h1_op=eps^(1/2)
h2=H(find(abs(df2-fp(x_0))==min(abs(df2-fp(x_0)))))
h2_op=(2*eps)^(1/(2+1))
h4=H(find(abs(df4-fp(x_0))==min(abs(df4-fp(x_0)))))
h4_op=(4*eps)^(1/(4+1))
h6=H(find(abs(df6-fp(x_0))==min(abs(df6-fp(x_0)))))
h6_op=(6*eps)^(1/(6+1))

%% Problem 4
close all
clear all
clc
% Part a
N = 1:100;
err = 0*N;
ff = @(x) abs(x).^3;
f = chebfun(ff,'splitting','on');
for k = 1:length(N)
    n = N(k);
    g = chebfun(ff,n);
    err(k) = norm(f-g,inf);
end
figure
loglog(N,err,'*',N,N.^-1,'--',N,N.^-3,'--')
set(gca,'fontsize',14)
grid on
xlabel('$N$ (log scale)','fontsize',20,'interpreter','latex')
ylabel('Error (log scale)','fontsize',20,'interpreter','latex')
saveas(gcf,'IMAGES/problem4a','epsc') 
% Part b
ff = @(x) exp(-1/sin(2*x^2));
f = chebfun(ff,'splitting','on');
for k = 1:length(N)
    n = N(k);
    g = chebfun(ff,n);
    err(k) = norm(f-g,inf);
end
figure
semilogy(N,err,'*',N,N.^-1,'--',N,N.^-3,'--')
set(gca,'fontsize',14)
grid on
xlabel('$N$','fontsize',20,'interpreter','latex')
ylabel('Error (log scale)','fontsize',20,'interpreter','latex')
saveas(gcf,'IMAGES/problem4b','epsc') 

% Part c
ff = @(x) 1/sin(1+x^2);
f = chebfun(ff,'splitting','on');
for k = 1:length(N)
    n = N(k);
    g = chebfun(ff,n);
    err(k) = norm(f-g,inf);
end
figure
semilogy(N,err,'*',N,N.^-1,'--',N,N.^-3,'--')
set(gca,'fontsize',14)
grid on
xlabel('$N$ ','fontsize',20,'interpreter','latex')
ylabel('Error (log scale)','fontsize',20,'interpreter','latex')
saveas(gcf,'IMAGES/problem4c','epsc') 

% Part d
ff = @(x) sinh(x)^2;
f = chebfun(ff,'splitting','on');
for k = 1:length(N)
    n = N(k);
    g = chebfun(ff,n);
    err(k) = norm(f-g,inf);
end
figure
semilogy(N,err,'*',N,N.^-1,'--',N,N.^-3,'--')
set(gca,'fontsize',14)
grid on
xlabel('$N$ (log scale)','fontsize',20,'interpreter','latex')
ylabel('Error (log scale)','fontsize',20,'interpreter','latex')
saveas(gcf,'IMAGES/problem4d','epsc')

%% Problem 5
close all
clear all
clc
ff = @(x) 1./(1+25*x.^2);
N = 1:100;
err = 0*N;
xx = linspace(-1,1,1000)';
for k = 1:length(N)
    x = linspace(-1,1,N(k))';
    W = baryWeights(x);
    y = ff(x);
    yy = bary(xx,y,x,W);
    err(k) = norm(ff(xx)-yy,inf);
end
figure
semilogy(N,err,'*',N,N.^-1,'--',N,N.^-3,'--')
set(gca,'fontsize',14)
grid on
xlabel('$N$','fontsize',20,'interpreter','latex')
ylabel('Error (log scale)','fontsize',20,'interpreter','latex')
saveas(gcf,'IMAGES/problem5_1','epsc') 

figure
plot(xx,ff(xx))
hold on
plot(xx,yy,'r--')
xlabel('$x$','fontsize',20,'interpreter','latex')
ylabel('$f(x),p(x)$','fontsize',20,'interpreter','latex')
legend({'$f(x)$','$p(x)$'},'Interpreter','latex','Location','north','fontsize',16)
axis([-1 1 -2 2])
grid on
saveas(gcf,'IMAGES/problem5_2','epsc') 