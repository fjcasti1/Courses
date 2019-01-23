%% HOMEWORK 2 - FRANCISCO CASTILLO

%%

%% Problem 2
clear all
close all
clc
format long
legendfontsize=12;
labelfontsize=14;
f = @(x) exp(sin(5*x));
f= chebfun(f);
f_cheb = chebfun(f,11);
xx=linspace(-1,1,1000);
P=legpoly(0:10,'norm');

figure
plot(P)
grid on
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
saveas(gcf,'Latex/FIGURES/LegendrePols','epsc')
saveas(gcf,'Latex/FIGURES/LegendrePols','fig')

f_N=0;
for k=1:11
    f_N=f_N+(f'*P(:,k))*P(:,k);
end
%%
figure
plot(f,'linewidth',2)
hold on
plot(f_N,'linewidth',2)
plot(f_cheb,'linewidth',2)
grid on
legend({'$f(x)=e^{\sin(5x)}$','Legendre','Chebishev'}...
    ,'fontsize',legendfontsize,'interpreter','latex','location','north')
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
% saveas(gcf,'Latex/FIGURES/Approximations_p2','epsc')
% saveas(gcf,'Latex/FIGURES/Approximations_p2','fig')

Error_cheb=norm(f-f_cheb,2)
Error_app=norm(f-f_N,2)
%%
figure
semilogy(abs(f-f_cheb))
hold on
semilogy(abs(f-f_N))
grid on
axis([-1 1 1e-4 0.35])
legend({'$e_{Cheb}$','$e_{Leg}$'},'fontsize',legendfontsize,'interpreter','latex')
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
% saveas(gcf,'Latex/FIGURES/Error_p2','epsc')
% saveas(gcf,'Latex/FIGURES/Error_p2','fig')

%% Problem 4
close all
N=10;
f = @(x) cos(3*x);
xx=linspace(-1,1,1000);
xk0=linspace(-1,1,N+1)';
xk1=chebpts(N+1,1);
xk2=chebpts(N+1,2);
l0 = @(x) prod(x-xk0);
l1 = @(x) prod(x-xk1);
l2 = @(x) prod(x-xk2);
%
figure
plot(xk0,l0(xk0),'r*')
hold on
plot(xx,l0(xx))
grid on
% legend({'Roots'},'fontsize',legendfontsize,'interpreter','latex')
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
saveas(gcf,'Latex/FIGURES/Lagrange','epsc')
saveas(gcf,'Latex/FIGURES/Lagrange','fig')

%
figure
plot(xk1,l1(xk1),'r*')
hold on
plot(xx,l1(xx))
grid on
% legend({'Roots'},'fontsize',legendfontsize,'interpreter','latex','location','north')
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
saveas(gcf,'Latex/FIGURES/Chebishev1','epsc')
saveas(gcf,'Latex/FIGURES/Chebishev1','fig')

%
figure
plot(xx,l2(xx))
hold on
plot(xk2,l2(xk2),'r*')
grid on
% legend({'Second kind Chebishev','Roots'},'fontsize',legendfontsize,'interpreter','latex')
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
saveas(gcf,'Latex/FIGURES/Chebishev2','epsc')
saveas(gcf,'Latex/FIGURES/Chebishev2','fig')

w0 = baryWeights(xk0);
p0 = @(x) bary(x,f(xk0),xk0,w0);
err_bound0 = 3^(N+1)*l0(xx)/factorial(N+1);
figure
semilogy(xx,abs(err_bound0))
hold on
semilogy(xx,abs(f(xx)-p0(xx)))
grid on
axis([-1 1 1e-10 1e-4])
legend({'$e_{Cauchy}$','$e_{equi}$'},'fontsize',legendfontsize,'interpreter','latex','Location','north')
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
ylabel('$e(x)$','fontsize',labelfontsize,'interpreter','latex')
saveas(gcf,'Latex/FIGURES/Cauchy_equi','epsc')
saveas(gcf,'Latex/FIGURES/Cauchy_equi','fig')

w1 = baryWeights(xk1);
p1 = @(x) bary(x,f(xk1),xk1,w1);
err_bound1 = 3^(N+1)*l1(xx)/factorial(N+1);
figure
semilogy(xx,abs(err_bound1))
hold on
semilogy(xx,abs(f(xx)-p1(xx)))
grid on
axis([-1 1 1e-10 1e-4])
legend({'$e_{Cauchy}$','$e_{Cheb1}$'},'fontsize',legendfontsize,'interpreter','latex','location','north')
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
ylabel('$e(x)$','fontsize',labelfontsize,'interpreter','latex')
saveas(gcf,'Latex/FIGURES/Cauchy_cheb1','epsc')
saveas(gcf,'Latex/FIGURES/Cauchy_cheb1','fig')

w2 = baryWeights(xk2);
p2 = @(x) bary(x,f(xk2),xk2,w2);
err_bound2 = 3^(N+1)*l2(xx)/factorial(N+1);
figure
semilogy(xx,abs(err_bound2))
hold on
semilogy(xx,abs(f(xx)-p2(xx)))
grid on
axis([-1 1 1e-10 1e-4])
legend({'$e_{Cauchy}$','$e_{Cheb2}$'},'fontsize',legendfontsize,'interpreter','latex','location','north')
xlabel('$x$','fontsize',labelfontsize,'interpreter','latex')
ylabel('$e(x)$','fontsize',labelfontsize,'interpreter','latex')
saveas(gcf,'Latex/FIGURES/Cauchy_cheb2','epsc')
saveas(gcf,'Latex/FIGURES/Cauchy_cheb2','fig')


E_equi=norm(f(xx)-p0(xx),2)
E_Cauchy_equi=norm(err_bound0,2)
E_Chebishev1=norm(f(xx)-p1(xx),2)
E_Cauchy_cheb1=norm(err_bound1,2)
E_Chebishev2=norm(f(xx)-p2(xx),2)
E_Cauchy_cheb2=norm(err_bound2,2)
