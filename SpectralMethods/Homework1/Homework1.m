clear variables; close all; clc
figformat='epsc';
%% Problem 4
N = 20:10:500;
x0 = -2*pi;
xf = 4*pi;
x = chebfun('x',[x0 xf]);
f = sin(x/2).^3;

figure
plot(f,'linewidth',2)
hold on
plot([0 0],[-1 1],'r-')
plot([2*pi 2*pi],[-1 1],'r-')
grid on
xlabel('$x$','interpreter','latex')
ylabel('$f(x)=\sin^3(x/2)$','interpreter','latex')
set(gca,'fontsize',14)
set(gca,'XTick',x0:pi/2:xf) 
xticklabels({'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi','7\pi/2','4\pi'})
txt='Latex/FIGURES/P4_1';
saveas(gcf,txt,figformat)

x0 = -2*pi;
xf = 0;
x = chebfun('x',[x0 xf]);
fp1 = -sin(x/2).^3;
x0 = 0;
xf = 2*pi;
x = chebfun('x',[x0 xf]);
fp2 = sin(x/2).^3;
x0 = 2*pi;
xf = 4*pi;
x = chebfun('x',[x0 xf]);
fp3 = -sin(x/2).^3;

x0 = -2*pi;
xf = 4*pi;
figure
plot(fp1,'b-','linewidth',2)
hold on
plot(fp2,'b-','linewidth',2)
plot(fp3,'b-','linewidth',2)
plot([0 0],[0 1],'r-')
plot([2*pi 2*pi],[0 1],'r-')
grid on
xlabel('$x$','interpreter','latex')
ylabel('$f_p(x)=\sin^3(x/2)$','interpreter','latex')
set(gca,'fontsize',14)
set(gca,'XTick',x0:pi/2:xf) 
xticklabels({'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi','7\pi/2','4\pi'})
txt='Latex/FIGURES/P4_2';
saveas(gcf,txt,figformat)
%%
N = 20:10:500;
x0 = 0;
xf = 2*pi;
x = chebfun('x',[x0 xf]);
f = sin(x/2).^3;

for j = 1:length(N)
    A = exp(1i*x*(-N(j):N(j)));
    lambda = 1/(2*pi)*A'*f;
    fn = A*lambda;
    err(j) = norm(f-fn,Inf);
end
%%
figure
loglog(err,'*','MarkerSize',12)
hold on
loglog(N,abs(N).^-3/1000,'r-')
grid on
xlabel('$N$','interpreter','latex')
ylabel('Error','interpreter','latex')
set(gca,'fontsize',14)
txt='Latex/FIGURES/P4_3';
saveas(gcf,txt,figformat)
ccc
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
