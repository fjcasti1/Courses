%% Problem 4
clear variables; close all; clc
figformat='epsc';

% Plot function
x0 = -2*pi;
xf = 2*pi;
x = chebfun('x',[x0 xf]);
f = sin(x/2).^3;

figure
plot(f,'b-','linewidth',2)
hold on
plot([0 0],[-1 1],'r--')
grid on
axis([x0 xf -1 1])
xlabel('$x$','interpreter','latex')
ylabel('$f(x)=\sin^3(x/2)$','interpreter','latex')
set(gca,'fontsize',14)
set(gca,'XTick',x0:pi/2:xf) 
xticklabels({'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'})
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

x0 = -2*pi;
xf = 2*pi;
figure
plot(fp1,'b-','linewidth',2)
hold on
plot(fp2,'b-','linewidth',2)
plot([0 0],[0 1],'r--')
grid on
axis([x0 xf 0 1])
xlabel('$x$','interpreter','latex')
ylabel('$f_p(x)=|\sin^3(x/2)|$','interpreter','latex')
set(gca,'fontsize',14)
set(gca,'XTick',x0:pi/2:xf) 
xticklabels({'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'})
txt='Latex/FIGURES/P4_2';
saveas(gcf,txt,figformat)

% Approximation and Error
N = 20:10:500;
x0 = 0;
xf = 2*pi;
x = chebfun('x',[x0 xf]);
f = sin(x/2).^3;

parfor j = 1:length(N)
    A = exp(1i*x*(-N(j):N(j)));
    lambda = 1/(2*pi)*A'*f;
    fn = A*lambda;
    err(j) = norm(f-fn,Inf);
end

% Plot Error
figure
loglog(N,err,'*','MarkerSize',12)
hold on
loglog(N,abs(N).^-3,'r-','linewidth',2)
grid on
% axis([20 510 1e-10 1e-1])
xlabel('$N$','interpreter','latex')
ylabel('Error','interpreter','latex')
legend('Approximation Error','Error Bound')
set(gca,'fontsize',14)
txt='Latex/FIGURES/P4_3';
saveas(gcf,txt,figformat)

%% Problem 5
clear variables; close all; clc
figformat='epsc';

% Plot function

x0 = -2*pi;
xf = 2*pi;
x = chebfun('x',[x0 xf]);
f = 3./(5-4*cos(x));

figure
plot(f,'b-','linewidth',2)
hold on
plot([0 0],[0 3.5],'r--')
grid on
xlabel('$x$','interpreter','latex')
ylabel('$f(x)=\frac{3}{5-4\cos(x)}$','interpreter','latex')
set(gca,'fontsize',14)
set(gca,'XTick',x0:pi/2:xf) 
xticklabels({'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'})
txt='Latex/FIGURES/P5_1';
saveas(gcf,txt,figformat)

%%
% Approximation and error
N = 20:10:500;
x0 = 0;
xf = 2*pi;
x = chebfun('x',[x0 xf]);
f = 3./(5-4*cos(x));

parfor j = 1:length(N)
    A = exp(1i*x*(-N(j):N(j)));
    lambda = 1/(2*pi)*A'*f;
    fn = A*lambda;
    err(j) = norm(f-fn,Inf);
end
%%
figure
loglog(N,err,'*','MarkerSize',12)
hold on
loglog(N,2.^(-abs(N)),'r','linewidth',2)
grid on
axis([10 500 1e-16 1e-1])
xlabel('$N$','interpreter','latex')
ylabel('Error','interpreter','latex')
legend('Approximation Error','Error Bound')
set(gca,'fontsize',14)
txt='Latex/FIGURES/P5_2';
saveas(gcf,txt,figformat)
