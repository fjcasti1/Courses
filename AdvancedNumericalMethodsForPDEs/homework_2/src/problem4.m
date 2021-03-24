close all; clear all; clc;

%% Global variables
global Npart L slow_interval v_slow v_fast
figure(1)
set(gcf, 'Position', get(0, 'Screensize'));
%% Parameter values
Npart = 500; % Number of particles
L = 50;
d = 0.0099; % Floating point error can cause the conditions to fail if exactly 0.01
T = 1;
dt = 1e-4;

% Fixed grid to compute density
M = Npart/2; % Number of grid cells
dy = L/M;
y = 0:dy:L;
u_max = (1/d+1/dy)/Npart; % Maximum posible density

% Speeds of the track
slow_interval = [30 40];
v_slow = 20;
v_fast = 100;

% Initial placement of particles
alpha= 0.02;
xi_0 = alpha*linspace(1,Npart,Npart);
% plot_particles(xi_0,[1:Npart],[])
plot_particles(xi_0,1:Npart,[])

% Prepare for the loop
xi = xi_0;
v = zeros(size(xi));
t = 0;
tstep = 0; % Time step counter

while t < T
    [idx_leaders,idx_followers] = find_leaders(xi, d, dt); % Array with the indeces of the leaders

    n = 1;
    for i=1:length(idx_leaders)
        idx = idx_leaders(i);
        v_leader = calculate_v0(mod(xi(idx),L));
        v(n:idx) = v_leader; % Match the speed of the leader
        n = idx + 1;
    end
    if (idx_leaders(end)<Npart) % This means the first particle is going to hit the last
        v(n:Npart) = v(1); % Match the velocities of the first and last particle
    end
    % Update xi, density u, t and tstep
    xi = mod(xi + dt*v,L);
    u  = calculate_density(y,dy,M,xi);

    t = t + dt;
	tstep = tstep + 1;
    
    % Plot solutions every 10 timesteps
    if (mod(tstep,20)==0)
        plot_particles(xi,idx_leaders,idx_followers)
        plot_density(y,u,t,u_max)
    end
end

function [leaders,followers] = find_leaders(xi, d, dt)
    global Npart L
    
    leaders = [];
	followers = [];
    
    v0 = calculate_v0(mod(xi,L)); % Free track velocity
    n = 1; % Particle counter
    
    % Index of the last one in the current lap
	last_in_lap = find(xi==max(xi)); 
    
    while n <= length(xi)
        % The particle in front of the last one is the first one
        if n == Npart
            next = 1;
        else
            next = n+1;
        end
        
        % Position if particle moved freely
        futurePosition = xi(n) + dt*v0(n);
        
        % Check if leader or follower and append to arrays
        if( n == last_in_lap && futurePosition + d <= xi(next) + L )
            leaders = [leaders n];
        elseif ( futurePosition + d <= xi(next) )
            leaders = [leaders n];
        else
            followers = [followers n];
        end
        
        % Update counter
        n = n+1;
    end
end

function v0 = calculate_v0(x)
    global slow_interval v_slow v_fast
    
    is_slow = x>slow_interval(1) & x<slow_interval(2);
    is_fast = ~is_slow;
    
    v0 = v_slow*is_slow + v_fast*is_fast;
end

function u = calculate_density(y,dy,M,xi)
    global Npart
    u = zeros(size(y));
    for m=1:M
        I=find((y(m)<xi)&(xi<y(m+1)));
        u(m)=length(I)/Npart/dy;
    end
end

function plot_particles(xi,idx_leaders,idx_followers)
    global L slow_interval
    
    % Plotting parameters
    linewidth = 2;
    
    % Radius of the circular track
    R = L/2*pi;
    
    % Coordinates of the track in the circle
    theta = linspace(0,2*pi,500);
    x = R*cos(theta);
    y = R*sin(theta);

    % Coordinates of first particle in the circle    
    theta_first = xi(1)*2*pi/L;
    x_first = R*cos(theta_first);
    y_first = R*sin(theta_first);
    % Coordinates of last particle in the circle    
    theta_last = xi(end)*2*pi/L;
    x_last = R*cos(theta_last);
    y_last = R*sin(theta_last);
    
    % Coordinates of leaders in the circle    
    theta_leaders = xi(idx_leaders(1:end))*2*pi/L;
    x_leaders = R*cos(theta_leaders);
    y_leaders = R*sin(theta_leaders);
    
    % Coordinates of followers in the circle
    theta_followers = xi(idx_followers(1:end))*2*pi/L;
    x_followers = R*cos(theta_followers);
    y_followers = R*sin(theta_followers);
    
    % Slow interval
    theta_slow = 2*pi*slow_interval/L;
    
    figure(1)
    subplot(1,2,1)
    % Plot track
    plot(x,y,'k', 'linewidth',linewidth/4)
    hold on
    
    % Plot followers
    plot(x_followers,y_followers,'r.')
    
    % Plot leaders
    plot(x_leaders,y_leaders,'b.')
    
    % Plot first and last particle
    plot(x_first,y_first,'gs','MarkerFaceColor','g')
    plot(x_last,y_last,'ms','MarkerFaceColor','m')
    
    % Plot slow area
    pl1 = line([1.05*R*cos(theta_slow(1)) 0],[1.05*R*sin(theta_slow(1)) 0],'linewidth',2);
    pl2 = line([1.05*R*cos(theta_slow(2)) 0],[1.05*R*sin(theta_slow(2)) 0],'linewidth',2);
    pl1.Color = 'red';
    pl1.LineStyle = '--';
    pl2.Color = 'red';
    pl2.LineStyle = '--';
    
    % Some options
    grid on
    pbaspect([1 1 1])
    hold off
end

function plot_density(y, u, t, u_max)
    global L slow_interval
    % Plotting parameters
    linewidth = 2;
    labelfontsize = 18;

    
    figure(1)
    subplot(1,2,2)
    plot(y,u, 'linewidth',linewidth)
    grid on
    xlabel('$x$','interpreter','latex','fontsize',labelfontsize)
    ylabel(['$u(x,t=',num2str(t),')$'],'interpreter','latex','fontsize',labelfontsize)
    axis([0 L 0 u_max])
    
    % Plot shaded area
    a = slow_interval(1);
    b = slow_interval(2);
    xx = [a a b b];
    yy = [0 u_max u_max 0];
    patch(xx,yy,'red')
    alpha(0.3)
    pbaspect([1 1 1])
end