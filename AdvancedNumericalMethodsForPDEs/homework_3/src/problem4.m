clear all, close all, clc

N = 50; % Number of terms in expansion
n = (-N:N)'; % (2N+1, 1)

% Define vector x
Nx = 100; % Number of grid steps
x = linspace(-pi,pi,Nx+1); % (Nx+1, 1)

% Construct matrix A
A = calculate_A(n,N);

% Calculate dt
dt = calculate_dt(A);

% Define useful matrices
M1 = eye(size(A)) - 0.5*dt*A; 
M2 = eye(size(A)) + 0.5*dt*A;
M = M1\M2; % As in c(k+1) = M c(k)

% Define initial coefficients for initial condition
c0 = calculate_initial_condition(n);
psi = calculate_wave_function(c0,n,x);
u = calculate_density(psi);
plot_psi(psi,x,0,0,N)
plot_u(u,x,0,0,N)

%% Find the time-steps to plot
T = 1;
t=0:dt:T;
t_plots = [0.33 0.66 0.99]; % Times to plot
tstep_plots = zeros(3,1);
for i=1:length(t_plots)
    [minValue, minIdx] = min(abs(t-t_plots(i)));
    tstep_plots(i) = minIdx;
end
tstep_plots
t(tstep_plots)
keyboard
%% Compute solutions and plot
for idx=1:length(tstep_plots)
    tstep = tstep_plots(idx);
    c = M^tstep*c0;
    t = (tstep-1)*dt;       
    psi = calculate_wave_function(c,n,x);
    u = calculate_density(psi);
    plot_psi(psi,x,t,idx,N)
    plot_u(u,x,t,idx,N)
end

function A = calculate_A(n,N)
    m = n;
    I = 1i * (1 - (-1).^(n'-m)) ./ (2 * pi * (n'-m));
    I(1:2*N+2:(2*N+1)^2) = 1/2; % Change diagonal values
    A = -1i * (I + diag(n.^2/2));
end

function dt = calculate_dt(A)
    spectral_radius = max(eig(A'*A));
    dt = 1/ sqrt(spectral_radius);
end

function c0 = calculate_initial_condition(n)
    c0 = 1i * exp(1i*n*pi/2) .* ( 1 - exp(1i*n*pi/2) ) ./ (sqrt(2*pi) * n);
    c0(n == 0) = sqrt(pi/8);
end

function psi = calculate_wave_function(c,n,x)
    psi = transpose(c) * exp(1i * n * x)/sqrt(2*pi);
end

function u = calculate_density(psi)
    u = abs(psi).^2;
end

function plot_psi(psi,x,t,idx,N)
    labelfontsize = 18;
    
    figure(1)
    plot(x,psi,'linewidth',2)
    xlim([-pi pi])
    xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
    ylabel(append('$\psi(x, t = ',num2str(t),')$'),'interpreter','latex','fontsize',labelfontsize)
    grid on
    figName = create_figName('psi',idx,N);
    saveas(gcf,figName,'png')
end

function plot_u(u,x,t,idx,N)
    labelfontsize = 18;
    
    figure(2)
    plot(x,u,'linewidth',2)
    xlim([-pi pi])
    xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
    ylabel(append('$u(x, t = ',num2str(t),')$'),'interpreter','latex','fontsize',labelfontsize)
    grid on
    figName = create_figName('u',idx,N);
    saveas(gcf,figName,'png')
end

function figName = create_figName(name,idx,N)
    path = '../figures/';
	figName = append(path,'p4_',name,'_N',num2str(N),'_sample',num2str(idx));
end