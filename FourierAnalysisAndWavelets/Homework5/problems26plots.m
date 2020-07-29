%% MATLAB CODE - FRANCISCO CASTILLO
%% Preliminary Commands
clear all
close all
clc
linewidth=1.6;
labelfontsize=18;
legendfontsize=12;
%% Problem 2
x=0:0.0001:1;
psi2=2.5*(heaviside(4*x)-2*heaviside(4*x-1)+heaviside(4*x-2));
psi21=-1*(heaviside(4*x-2)-2*heaviside(4*x-3)+heaviside(4*x-4));
psi=-(5/4)*(heaviside(2*x)-2*heaviside(2*x-1)+heaviside(2*x-2));
phi=(3/4)*(heaviside(x));

figure
plot(x,psi2,'linewidth',linewidth);
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
title('$\frac{5}{2}\psi\left(2x\right)$','interpreter','latex','fontsize',labelfontsize)
axis([0 1 -4 4])
grid on
saveas(gcf,'psi2p2','png')

figure
plot(x,psi21,'linewidth',linewidth);
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
title('$-\psi\left(2x-1\right)$','interpreter','latex','fontsize',labelfontsize)
axis([0 1 -4 4])
grid on
saveas(gcf,'psi21p2','png')

figure
plot(x,psi,'linewidth',linewidth);
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
title('$-\frac{5}{4}\psi\left(x\right)$','interpreter','latex','fontsize',labelfontsize)
axis([0 1 -4 4])
grid on
saveas(gcf,'psip2','png')

figure
plot(x,phi,'linewidth',linewidth);
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
title('$\frac{3}{4}\phi\left(x\right)$','interpreter','latex','fontsize',labelfontsize)
axis([0 1 -4 4])
grid on
saveas(gcf,'phip2','png')

figure
plot(x,psi2+psi21+psi+phi,'linewidth',linewidth);
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
title('$f(x)=\frac{5}{2}\psi\left(2x\right)-\psi\left(2x-1\right)-\frac{5}{4}\psi\left(x\right)+\frac{3}{4}\phi\left(x\right)$','interpreter','latex','fontsize',labelfontsize)
axis([0 1 -4 4])
grid on
saveas(gcf,'fp2','png')
%% Problem 6
x=0:0.0001:1;
phi=0;
c=[-1 2 1 3 3 2 -2 -1]
for j=0:7
    phi=phi+c(j+1)*(heaviside(8*x-j)-heaviside(8*x-j-1));
end

figure
plot(x,phi,'linewidth',linewidth);
xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
ylabel('$g(x)$','interpreter','latex','fontsize',labelfontsize)
axis([0 1 -4 4])
grid on
saveas(gcf,'gp6','png')
%%
% h1=legend('$m=0$','$m=6$','$m=10$','$m=15$');
% set(h1,'interpreter','latex','fontsize',legendfontsize);
