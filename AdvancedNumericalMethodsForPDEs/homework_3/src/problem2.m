clear all, close all, clc


% Number of intervals
M = 100; % In x-direction
N = 100; % In Y-direction

% Define limits of domain
x_0 = 0;
x_F = 1;
y_0 = 0;
y_F = 1;
T = 0.08;

% Step sizes
dx = (x_F - x_0) / M;
dy = (y_F - y_0) / N;
dt = 1e-4;

% Define grid (including ghost nodes)
x = x_0 - dx : dx : x_F + dx; % x(1) and x(M+3) are ghost nodes
y = y_0 - dy : dy : y_F + dy; % y(1) and y(M+3) are ghost nodes

% Define useful constants
C = dt/dx;
mu_x = dt/dx^2;
mu_y = dt/dy^2;

% Define initial condition
U0 = kron((x==0.5),(y==0.5))';
plot_sol(U0,x,y,0)

% Define big matrices size, including ghost nodes
S = (M + 3)*(N + 3);

%% Construct Matrix A
% Diagonal component
d = diag((1 + mu_x + mu_y) * ones(S,1));
% 1-lower diagonal component
l1 = diag(-(C/4 + mu_x/2) * ones(S-1,1), -1);
% 1-upper diagonal component
u1 = diag((C/4 - mu_x/2) * ones(S-1,1), 1);
% (M+3)-lower diagonal component
lM3 = diag(-mu_x/2 * ones(S-M-3,1), -M-3);
% (M+3)-upper diagonal component
uM3 = diag(-mu_y/2 * ones(S-M-3,1), M+3);
% Add them together
A = sparse(d + u1 + l1 + lM3 + uM3);

%% Construct Matrix B
% Diagonal component
d = diag((1 - mu_x - mu_y) * ones(S,1));
% 1-lower diagonal component
l1 = diag((C/4 + mu_x/2) * ones(S-1,1), -1);
% 1-upper diagonal component
u1 = diag((mu_x/2 - C/4) * ones(S-1,1), 1);
% (M+3)-lower diagonal component
lM3 = diag(mu_x/2 * ones(S-M-3,1), -M-3);
% (M+3)-upper diagonal component
uM3 = diag(mu_y/2 * ones(S-M-3,1), M+3);
% Add them together
B = sparse(d + u1 + l1 + lM3 + uM3);

%% Compute solution

t_plots = [0.01 0.02 0.04 0.06 0.08];

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
%     U(i + 2*(M+3)) = 1;
    U(i) = U(i + 2*(M+3)); % Bottom
%     U(i + N*(M+3)) = 2;
    U(i + (N+2)*(M+3)) = U(i + N*(M+3)); % Top
    
%     U(3 + (M+3)*(j-1)) = 3/2
%     U(2 + (M+3)*(j-1)) = -3/(2*dx)/2
    U(1 + (M+3)*(j-1)) = U(3 + (M+3)*(j-1)) - 2*dx*U(2 + (M+3)*(j-1)); % Left
    
%     U(M+1 + (M+3)*(j-1)) = 4/2
%     U(M+2 + (M+3)*(j-1)) = 4/(2*dx)/2
    U(M+3 + (M+3)*(j-1)) = U(M+1 + (M+3)*(j-1)) + 2*dx*U(M+2 + (M+3)*(j-1)); % Right
    % Advance time
    t = round(t + dt,4); % avoid floating point error to compare with ismember
    tstep = tstep + 1;
    
    if (ismember(t,t_plots))
        plot_sol(U,x,y,t)
        pause(0.5)
    end
    
end


function plot_sol(U,x,y,t)
    
    Lx = length(x);
    Ly = length(y);
    u = reshape(U, [Lx Ly]);
    figure
    
    surf(x(2:Lx-1),y(2:Ly-1),u(2:Lx-1,2:Ly-1))
%     zlim([-0.5 1])
    ylim([0 1])
    xlim([0 1])
    xlabel('x')
    ylabel('y')
    zlabel('u(x,y,t)')
    title(num2str(t))
end