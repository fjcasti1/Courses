clear variables; close all; clc
figformat='epsc';
%% Problem 4
N = 20:10:500;
x = chebfun('x',[0 2*pi]);
f = sin(x/2).^3;

figure
plot(f,'linewidth',2)
grid on
xlabel('$x$','interpreter','latex')
ylabel('$f(x)=\sin^3(x/2)$','interpreter','latex')
set(gca,'fontsize',14)
txt='Latex/FIGURES/P4_1';
saveas(gcf,txt,figformat)

for j = 1:length(N)
    A = exp(1i*x*(-N(j):N(j)));
    lambda = 1/(2*pi)*A'*f;
    fn = A*lambda;
    err(j) = norm(f-fn,Inf);
end

figure
loglog(err,'*','MarkerSize',12)
grid on
xlabel('$N$','interpreter','latex')
ylabel('Error','interpreter','latex')
set(gca,'fontsize',14)
txt='Latex/FIGURES/P4_2';
saveas(gcf,txt,figformat)

%% Problem 5
N = 20:10:500;
x = chebfun('x',[0 2*pi]);
f = 3./(5-4*cos(x));

figure
plot(f,'linewidth',2)
grid on
xlabel('$x$','interpreter','latex')
ylabel('$f(x)=\frac{3}{5-4\cos(x)}$','interpreter','latex')
set(gca,'fontsize',14)
txt='Latex/FIGURES/P5_1';
saveas(gcf,txt,figformat)

for j = 1:length(N)
    A = exp(1i*x*(-N(j):N(j)));
    lambda = 1/(2*pi)*A'*f;
    fn = A*lambda;
    err(j) = norm(f-fn,Inf);
end

figure
loglog(err,'*','MarkerSize',12)
grid on
xlabel('$N$','interpreter','latex')
ylabel('Error','interpreter','latex')
set(gca,'fontsize',14)
txt='Latex/FIGURES/P5_2';
saveas(gcf,txt,figformat)
