clear all; close all; clc
format long

enableVideo = false;

a = 0.5;
b = 0;
dx = 0.1;
dt_default = calculate_dt(a,b,dx);
dt = dt_default;

x0 = -1;
xN = 1;
N = (xN-x0)/dx;
x = linspace(x0,xN,N+1)';
u0 = heaviside(-x);

[A,B,C] = calculateDiagonals(a,b,dx,dt);
Mtilde_default = calculate_Mtilde(A,B,C,N);
Mtilde = Mtilde_default;

t = 0;
u = u0;
T = 4;

k = 1:5;
plotTimes = k*T/5

storeCounter = 1;
shouldStore = false;
storedSolutions = [];
dtHasChanged = false;

while t<T
    if (dtHasChanged) % Need to reset values
        dt = dt_default;
        [A,B,C] = calculateDiagonals(a,b,dx,dt);
        Mtilde = Mtilde_default;
        dtHasChanged = false;
    end
    if(t+dt > plotTimes(storeCounter))
        dt = plotTimes(storeCounter) - t;
        % We need to recalculate the matrix for the new dt
        [A,B,C] = calculateDiagonals(a,b,dx,dt);
        Mtilde = calculate_Mtilde(A,B,C,N);
        shouldStore = true; % Should plot the solution after this iteration  
        dtHasChanged = true;
    end
    
    % -- Advance solution --
    u_prev = u;
    u(2:N) = Mtilde * u_prev; % Solve the interior
    % Periodic BCs
    u(1)   = A*u_prev(N) + B*u_prev(1) + C*u_prev(2);
    u(N+1) = A*u_prev(N) + B*u_prev(N+1) + C*u_prev(2);
    
    % -- Advance time -- 
    t = t + dt;
    
    if(shouldStore)
        disp(['Storing solution at t = ',num2str(t)])
        storedSolutions = [storedSolutions u];
        storeCounter = storeCounter +1;
        shouldStore = false;
    end
    if(enableVideo)
        figure(1)
        grid on
        plot(x,u);
        axis([-1 1 min(u0) max(u0)])
    end
end

figName = create_figName(b,dx);
plot_solutions(x,[0 plotTimes],[u0 storedSolutions],[-1 1 min(u0) max(u0)],figName)

function figName = create_figName(b,dx)
    figName = 'sol_b';
    if (b==0)
        figName = append(figName,'0_dx');
    else
        exponent = floor(log10(b));
        base = b/10^(exponent);
        figName = append(figName,num2str(base),'e',num2str(exponent),'_dx');
    end
    exponent = floor(log10(dx));
    base = dx/10^(exponent);
    figName = append(figName,num2str(base),'e',num2str(exponent)); 
end

function Mtilde = calculate_Mtilde(A,B,C,N)
    M = diag(A*ones(1,N),-1) + diag(B*ones(1,N+1)) + diag(C*ones(1,N),1);
    Mtilde = M(2:N,:); % For the interior
end

function [A,B,C] = calculateDiagonals(a,b,dx,dt)
    c = dt/dx; % Courant Number
    if (b==0) % Lax-Friedrichs
        A = (1 + a*c)/2; % Lower diagonal
        B = 0;           % Main diagonal
        C = (1 - a*c)/2; % Upper diagonal
    else % FTCS
        A = c*(b/dx + a/2); % Lower diagonal
        B = 1 - 2*b*c/dx;   % Main diagonal
        C = c*(b/dx - a/2); % Upper diagonal
    end
end

function plot_solutions(x,times,solutions,axisLimits,figName)
    linewidth = 2;
    labelfontsize = 18;
    legendfontsize = 12;
    
    figure(2)
    grid on
    hold on
    for i=1:length(times)
        plot(x,solutions(:,i),'DisplayName',['t = ',num2str(times(i))],'linewidth',linewidth);        
    end
    xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
    ylabel('$u(x,t)$','interpreter','latex','fontsize',labelfontsize)
    l = legend;
    set(l,'fontsize',legendfontsize)
    axis(axisLimits)
    saveas(gcf,figName,'png')
end

function round_number = round_down(number, decimals)
    multiplier = 10^decimals;
    round_number = floor(number * multiplier)/multiplier;
end

function dt = calculate_dt(a,b,dx)
    if (b==0) % Lax-Friedrichs
        dt = round_down(dx/abs(a), 6);
    else % FTCS
        dt = round_down(min(2*b/a^2, dx^2/(2*b)), 6);        
    end    
end