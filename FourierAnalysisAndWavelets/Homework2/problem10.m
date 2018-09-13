%% Problem 10
clear all
close all
clc

linewidth = 2;
fontsize = 15;
N=500;
x=linspace(0,pi,N);
x_ext=linspace(-pi,pi,2*N);
f = pi-x;
f_even=pi-x_ext;

figure
plot(x_ext(N+1:2*N),f,x_ext(1:N),-fliplr(f),'Linewidth',linewidth)
xlabel('$x$','Interpreter','latex','Fontsize',fontsize)
ylabel('$f(x)$','Interpreter','latex','Fontsize',fontsize)
axis([-pi pi -pi pi])
grid on
legend('Original function','Odd Extension')
saveas(gcf,'oddExtension','png')

figure
plot(x_ext(N+1:2*N),f,x_ext(1:N),fliplr(f),'Linewidth',linewidth)
xlabel('$x$','Interpreter','latex','Fontsize',fontsize)
ylabel('$f(x)$','Interpreter','latex','Fontsize',fontsize)
axis([-pi pi 0 pi])
grid on
legend('Original function','Even Extension')
saveas(gcf,'evenExtension','png')

figure
plot(x_ext(N+1:2*N),f,x_ext(1:N),-fliplr(f)+pi,'Linewidth',linewidth)
xlabel('$x$','Interpreter','latex','Fontsize',fontsize)
ylabel('$f(x)$','Interpreter','latex','Fontsize',fontsize)
axis([-pi pi 0 pi])
grid on
h=legend('Original function','$\pi$ Extension');
set (h,'Interpreter','latex')
saveas(gcf,'piExtension','png')