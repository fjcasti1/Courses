clear all; close all; clc
format long

enableVideo = false;
%% Parameter Values
L = 50;
c = 100;
N = 100; % Mesh size
T = 1;   % Final time
x0 = 0;
xN = L;
dx = (xN-x0)/N;
x = linspace(x0,xN,N+1)';

%% Initial condition and Free velocity
% Initial Condition
step_width = 10;
initial_speed = 50;
u0 = calculate_u0(x,initial_speed,step_width);
% Free velocity v0 (density u=0)
v_slow_interval = [30, 40];
v_slow = 20;
v_fast = 100;
v0 = calculate_v0(x,v_slow_interval,v_slow,v_fast);


%% Prepare for the loop
t = 0;
u = u0;
u_new = zeros(size(u)); % Just to allocate for its size

saveTimes = 0.2:0.2:T;
save_idx = 1;

% Plot initial condition
plot_solution(x,t,u,c,L,v_slow_interval)
saveas(gcf,'../figures/u_t0','png')

while t<T
    v_g = calculate_group_speed(v0,u,c); % Speed of the particles
    flux = u.*v_g;
    dt = calculate_dt(v_g,dx);
    
    % Check for plotting times
    if ( saveTimes & t+dt > saveTimes(save_idx)) %#ok<AND2,BDSCI,BDLGI>
        dt = saveTimes(save_idx) -t;
        save_idx = save_idx +1;
    elseif ( t+dt > T)
        dt = T -t;
    end
    
    C = dt/dx; % Courant Number
    % Advance in time
    for j=2:N % just the interior
        u_new(j) = ( u(j+1) + u(j-1) )/2 - C * ( flux(j+1) - flux(j-1) )/2;
    end
    
    % Periodict Boundary Conditions
    u_new(1) = ( u(2) + u(N+1) )/2 - C * ( flux(2) - flux(N+1) )/2; % BC on the left
    u_new(N+1) = ( u(1) + u(N) )/2 - C * ( flux(1) - flux(N) )/2; % BC on the right
    
    % Update u and t  
    u = u_new; 
    t = t +dt;
    
    % Plot solution and save if appropriate
    plot_solution(x,t,u,c,L,v_slow_interval)
    if (ismember(t,saveTimes))
        figName = create_figName(t);
        saveas(gcf,figName,'png')
        disp(['Saved figure -> "',figName,'.png"'])
    end
end

function v = calculate_group_speed(v0,u,c) % Speed of the particles
    v = v0.*(1-u/c);
end

function u0 = calculate_u0(x,initial_speed,step_width) % Initial solution
    u0 = initial_speed * heaviside(step_width - x);
end

function v0 = calculate_v0(x,slow_interval,v_slow,v_fast) % Free velocity (density u=0)
    is_slow = x>slow_interval(1) & x<slow_interval(2);
    is_fast = ~is_slow;
    v0 = v_slow*is_slow + v_fast*is_fast;
end

function dt = calculate_dt(v,dx)
    a = max(v);
    dt = 0.9*round_down(dx/abs(a),6);
end

function round_number = round_down(number,decimals)
    multiplier = 10^decimals;
    round_number = floor(number * multiplier)/multiplier;
end

function plot_solution(x,t,u,c,L,v_slow_interval)
    % Plot parameters
   a = v_slow_interval(1);
   b = v_slow_interval(2);
   linewidth = 2;
   labelfontsize = 18;

   figure(1)
   plot(x,u,'linewidth',linewidth)
   grid on
   axis([0 L 0 c])
   xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
   ylabel(['$u(x,t=',num2str(t),')$'],'interpreter','latex','fontsize',labelfontsize)

   % Plot shaded area
   xx = [a a b b];
   yy = [0 c c 0];
   patch(xx,yy,'red')
   alpha(0.3)
end

function figName = create_figName(t)
    exponent = floor(log10(t));
    base = t/10^exponent;
    path = '../figures/';
    figName = append(path,'u_t',num2str(base),'e',num2str(exponent));
end