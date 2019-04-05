%% Homework 4, Problem 1 - Francisco Castillo'
clear all; close all; clc;
labelfontsize = 14;

NN = 4:2:50;
% for j=1:length(NN)
    N=NN(j);
    [D,x] = cheb(N); % D:(N+1)x(N+1), x:(N+1)x1
    D2 = D^2;
    D2 = D2(2:N,2:N);       % Zero Dirichlet BCs
    E = diag(exp(u(2:N)));
    u = (D2+4*D+E)\f;
    u = [0;u;0];
    u0(j) = u(N/2+1);
% end

figure
plot(NN,u0,'r*')
grid on
xlabel('$N$','fontsize',labelfontsize,'interpreter','latex')
ylabel('$u(0)$','fontsize',labelfontsize,'interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
txt = 'Latex/FIGURES/P1_a';
saveas(gcf,txt,'png')

xx = [x(end):.01:x(1)]';
uu = polyval(polyfit(x,u,N),xx);

figure
plot(x,u,'b*')
hold on
plot(xx,uu,'b')
grid on
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
ylabel('$u(x)$','fontsize',labelfontsize,'interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
txt = 'Latex/FIGURES/P1_b';
saveas(gcf,txt,'png')
