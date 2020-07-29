%% HOMEWORK 2 - FRANCISCO CASTILLO
%
clear all; close all; clc;
%% Problem 2 
x0=0;
xf=2*pi;
N=11;
K=5;
L=1e3;
x=linspace(x0,xf,N+1);
x(end)=[];
k=-K:K;
l=-L:L;

fhat = 1./((1+1i)*(k+(1/2)).^6);
fxj = 1./((1+1i)*(l+(1/2)).^6)*exp(1i*l'*x);
vhat = (1/N)*fxj*exp(-1i*x'*k);

lhs = vhat-fhat;
m = l;
m(1001) = [];

rhs = zeros(1,11);
for l = 1:length(k)
   rhs(l) = sum(1./((1+1i)*((k(l)+m*N)+(1/2)).^6));
end
diff = lhs-rhs
%% Problem 3
tmax = 1;
eps = [1e-1 1e-2 1e-3]';
tol = 1e-3;
mesh = zeros(length(eps),1);
E = zeros(length(eps),1);
for m=1:length(eps)
    j = 5;
    N = 2^j;
    error = 1;
    u = PDE_solve(N,tmax,eps(m),false);
    while (error>tol && N<2048)
        j = j+1;
        N = 2^j;
        ufine = PDE_solve(N,tmax,eps(m),false);
        error = norm(u(end,:)-ufine(end,2:2:end),inf);
        u = ufine;
        j
    end
    m
    mesh(m) = N;
    E(m) = error;
end
mesh
E

%% Problem 4
close all
figformat='epsc';
linewidth=2;
% Compute derivatives for various values of N:
Nmax = 50; E = zeros(3,Nmax/2-2);
for N = 6:2:Nmax
    h = 2*pi/N; x = h*(1:N)';
    column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
    D = toeplitz(column,column([1 N:-1:2]));
    v = abs(sin(x)).^3;                     % 3rd deriv in BV
    vprime = 3*sin(x).*cos(x).*abs(sin(x));
    E(1,N/2-2) = norm(D*v-vprime,inf);
    v = exp(-sin(x/2).^(-2));               % C-infinity
    vprime = .5*v.*sin(x)./sin(x/2).^4;
    E(2,N/2-2) = norm(D*v-vprime,inf);
    v = 1./(1+sin(x/2).^2);                 % analytic in a strip
    vprime = -sin(x/2).*cos(x/2).*v.^2;
    E(3,N/2-2) = norm(D*v-vprime,inf);
    v = sin(10*x); vprime = 10*cos(10*x);   % band-limited
    E(4,N/2-2) = norm(D*v-vprime,inf);
end
Nvector = 6:2:Nmax;
hvector = 2*pi./Nvector;


p = 3;
v = 1;

figure
semilogy(Nvector,E(1,:),'*')
hold on
semilogy(Nvector,hvector.^(p-v),'linewidth',linewidth)
grid on
xlabel('$N$','interpreter','latex')
ylabel('$|w_j-u^{(v)(x_j)}|$','interpreter','latex')
set(gca,'fontsize',14)
txt='Latex/FIGURES/P4_1';
saveas(gcf,txt,figformat)

a = 1.76;
e = 1e-4;

figure
semilogy(Nvector,E(3,:),'*')
hold on
semilogy(Nvector,exp(-pi*(a-e)./hvector),'linewidth',linewidth)
grid on
xlabel('$N$','interpreter','latex')
ylabel('$|w_j-u^{(v)(x_j)}|$','interpreter','latex')
set(gca,'fontsize',14)
txt='Latex/FIGURES/P4_2';
saveas(gcf,txt,figformat)

figure
semilogy(Nvector,E(4,:),'*')
hold on
plot([20 20],[1e-16 10],'r--','linewidth',linewidth)
axis([5 51 1e-16 2e1])
grid on
xlabel('$N$','interpreter','latex')
ylabel('$|w_j-u^{(v)(x_j)}|$','interpreter','latex')
set(gca,'fontsize',14)
txt='Latex/FIGURES/P4_3';
saveas(gcf,txt,figformat)