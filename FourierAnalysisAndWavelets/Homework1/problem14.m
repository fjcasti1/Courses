clear all; close all; clc
format long
x=linspace(-pi,pi,200);
f0=pi^2/3-4*cos(x)+cos(2*x);
g0=2*(pi^2-6)*sin(x)+(3/2-pi^2)*sin(2*x);

plot(x,f0,x,g0,'linewidth',1.5)
axis([-pi pi -25 25])
xlabel('$x$','Interpreter','latex')
grid on
legend('f_0','g_0','Interpreter','latex')
saveas(gcf,'projectionsP14','epsc')