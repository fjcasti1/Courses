clear all, close all, clc

format short

% Number of intervals
M = 20; % In x-direction
N = 20; % In Y-direction

% Define limits of domain
x_0 = 0;
x_F = 1;
y_0 = 0;
y_F = 1;
T = 0.2;

% Step sizes
dx = (x_F - x_0) / M;
dy = (y_F - y_0) / N;
dt = 0.0001;

% Define grid (including ghost nodes)
x = linspace (x_0-dx, x_F+dx, M+3); % x(1) and x(M+3) are ghost nodes
y = linspace (y_0-dy, y_F+dy, N+3); % y(1) and y(N+3) are ghost nodes

% Define useful constants
C = dt/dx;
mu_x = dt/(dx^2);
mu_y = dt/(dy^2);

% Define initial condition
U0 = kron((round(x,8) == 0.5),((round(y,8) == 0.5)))'/dx/dy;
plot_sol(U0,x,y,0)

% Define big matrices size, including ghost nodes
S = (M + 3)*(N + 3);

%% Construct Matrix A
% Diagonal component
d = diag((1 + mu_x + mu_y) * ones(S,1));
% 1-lower diagonal component
l1 = diag((-C/4 - mu_x/2) * ones(S-1,1), -1);
% 1-upper diagonal component
u1 = diag(( C/4 - mu_x/2) * ones(S-1,1), 1);
% (M+3)-lower diagonal component
lM3 = diag(-mu_y/2 * ones(S-M-3,1), -M-3);
% (M+3)-upper diagonal component
uM3 = diag(-mu_y/2 * ones(S-M-3,1), M+3);
% Add them together
A = sparse(d + u1 + l1 + lM3 + uM3);

%% Construct Matrix B
% Diagonal component
d = diag((1 - mu_x - mu_y) * ones(S,1));
% 1-lower diagonal component
l1 = diag(( C/4 + mu_x/2) * ones(S-1,1), -1);
% 1-upper diagonal component
u1 = diag((-C/4 + mu_x/2) * ones(S-1,1), 1);
% (M+3)-lower diagonal component
lM3 = diag(mu_y/2 * ones(S-M-3,1), -M-3);
% (M+3)-upper diagonal component
uM3 = diag(mu_y/2 * ones(S-M-3,1), M+3);
% Add them together
B = sparse(d + u1 + l1 + lM3 + uM3);

%% Compute solution

t_plots = 0.01:0.01:T

t = 0;
tstep = 0;
plot_idx = 1;
% Initialize
U = U0;

i = 1:M+3; % Horizontal array index
j = 1:N+3; % Vertical array index
while t < T
        
    % Advance solution
    rhs = B*U;
    U = A\rhs;
    
    % Impose boundary conditions
    U(i) = U(i + 2*(M+3)); % Bottom
    U(i + (N+2)*(M+3)) = U(i + N*(M+3)); % Top
    U(1 + (M+3)*(j-1)) = U(3 + (M+3)*(j-1)) - 2*dx*U(2 + (M+3)*(j-1)); % Left
    U(M+3 + (M+3)*(j-1)) = U(M+1 + (M+3)*(j-1)) + 2*dx*U(M+2 + (M+3)*(j-1)); % Right
    
    % Advance time
    t = round(t + dt,4); % avoid floating point error to compare with ismember
    tstep = tstep + 1;
    
    if (mod(tstep,100) == 0)
        plot_sol(U,x,y,t)
        pause(0.15)
    end
    
end

function plot_sol(U,x,y,t)
    labelfontsize = 18;
    
    Lx = length(x);
    Ly = length(y);
    u = reshape(U, [Lx Ly]);
    figure(1)
    
    surf(x(2:Ly-1),y(2:Lx-1),u(2:Lx-1,2:Ly-1)')

    xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
    ylabel('$y$','interpreter','latex','fontsize',labelfontsize)
    zlabel(append('$u(x, y, t = ',num2str(t),')$'),'interpreter','latex','fontsize',labelfontsize)
    figName = create_figName(t);
    saveas(gcf,figName,'png')
end

function figName = create_figName(t)
    if t == 0
        exponent = 0;
        base = 00;
    else
        exponent = -2;
        base = t*100;
    end
    
    path = '../figures/';
    
    if base > 9
        figName = append(path,'p2_u_t',num2str(base),'e',num2str(exponent));
    else
        figName = append(path,'p2_u_t0',num2str(base),'e',num2str(exponent));
    end
end