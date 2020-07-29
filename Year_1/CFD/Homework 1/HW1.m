%% HOMEWORK 1 - FRANCISCO CASTILLO

%% Problem 2
format long;clear all;close all;clc
x_0=33.3;
f = @(x) (1+2*x^2*cos(x))/(x^(2.4));
fp = @(x) (2*x^2*(2*cos(x)-x*sin(x))-2.4*(1+2*x^2*cos(x)))/(x^3.4);
for k = -3:1:25
    H(k+4) = 2^(-k);
end
for k = 1:length(H)
    h = H(k);
    df1(k) = (f(x_0+h)-f(x_0))/h;
    df2(k) = (f(x_0+h)-f(x_0-h))/(2*h);
    df6(k) = (45*(f(x_0+h)-f(x_0-h))-9*(f(x_0+2*h)-f(x_0-2*h))...
        +(f(x_0+3*h)-f(x_0-3*h)))/(60*h);
end
linewidth=2;
darkgreen=[0 0.6 0];
figure('units','normalized','outerposition',[0 0 1 1])
loglog(H,abs(df1-fp(x_0)),'*',H,abs(df2-fp(x_0)),'*')
hold on
loglog(H,abs(df6-fp(x_0)),'*','Color',darkgreen)
loglog(H,H/2,'b--',H,H.^(2)/6,'r--','linewidth',linewidth)
loglog(H,H.^(6)/140,'--','Color',darkgreen,'linewidth',linewidth)
loglog(H,eps./H,'k:','linewidth',linewidth)
set(gca,'fontsize',14)
axis([1e-8 1e1 1e-20 1e5])
grid on
xlabel('$h$ (log scale)','fontsize',20,'interpreter','latex')
ylabel('Error (log scale)','fontsize',20,'interpreter','latex')
saveas(gcf,'IMAGES/problem2_1','epsc')
saveas(gcf,'IMAGES/problem2_1','fig')