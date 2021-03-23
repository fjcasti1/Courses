close all; clear all; clc;
format long
%% Plotting parameters
global linewidth
linewidth = 2;

%% Parameter values
global Npart L d slow_interval

Npart = 400;
L = 50;
d = 0.0099;
slow_interval = [39,40];
T = 200;
dt = 1e-4;
alpha = L/Npart;
% alpha= 0.02;

% Initial placement of particles
% xi_0 = 0.9*alpha*linspace(1,Npart,Npart);
% xi_0 = sort(rand(1,Npart)*L);
xi_0 = alpha*linspace(1,Npart,Npart);
% xi_0(end) = 50.00995 - 0.09996;
% xi_0(end-1) = xi_0(end)-0.1;
% xi_0(end-2) = xi_0(end-1)-0.1;
plot_particles(xi_0,[],[])


% Prepare for the loop
xi = xi_0;

v = zeros(size(xi));
t = 0;
tstep = 0;
n = 1;

while t < T
	
    [idx_leaders,idx_followers] = find_leaders(xi, dt); % Array with the indeces of the leaders
    n = 1;
    for i=1:length(idx_leaders)
        idx = idx_leaders(i);
        v(n:idx) = calculate_v0(mod(xi(idx),L));
        n = idx + 1;
    end
    if (idx_leaders(end)<Npart) % This means the first particle is going to hit the last
        v(n:Npart) = v(1); % Match the velocities of the first and last particle
    end
%     if(idx_leaders(end)<Npart)
%         disp('last leader')
%         idx_leaders(end)
%     end
    xi = mod(xi + dt*v,L);
    if (xi(end)+5*dt*v(end)>xi(1)+L)
        keyboard
    end
    t = t + dt;
%     xi(end)
%     xi(1)
%     [M,I] = max(xi)
    if (mod(tstep,10)==0)
        plot_particles(xi,idx_leaders,idx_followers)
    end
%     idx_leaders
    tstep = tstep + 1;
end

function [leaders,followers] = find_leaders(xi, dt)
    global Npart L d ;
    
    leaders = [];
	followers = [];
    
    n = 1;
    
    while n <= length(xi)
        v0 = calculate_v0(mod(xi(n),L));
        
        % The particle in front of the last one is the first one
        if n == Npart
            next = 1;
        else
            next = n+1;
        end
        
        futurePosition = xi(n) + dt*v0;
%         if (n>292)
%             keyboard
%         end
     
%         xi
%         max(xi)
        last_in_lap = n == find(xi==max(xi)); % Index of the last one in the current lap
%         if last_in_lap
%             n
%             last_in_lap
%         end
%         fprintf('Future Position = %2.6f\n',futurePosition)
%         fprintf('Next Particle = %2.6f\n',xi(next))

        % Check if we are too close to the next particle to calculate velocity
%         if futurePosition>L
%         
%             keyboard
%         end
        if( last_in_lap && futurePosition + d <= xi(next) + L )
            leaders = [leaders n];
        elseif ( futurePosition + d <= xi(next) )
            leaders = [leaders n];
        else
            followers = [followers n];
        end
        
        n = n+1;
        
    end
%     R = L/2*pi;
%     for n=1:length(leaders)    
%         theta_part = xi(leaders(n))*2*pi/L;
%         x_part = R*cos(theta_part);
%         y_part = R*sin(theta_part); 
%         
%         figure(1)
%         hold on
%         plot(x_part,y_part,'b*')
%         grid on
%         pbaspect([1 1 1])
%     end
end

function v0 = calculate_v0(x)
    global slow_interval
    if (x>slow_interval(1) && x<slow_interval(2))
        v0 = 10;
    else
        v0 = 100;
    end
end

function plot_particles(xi,idx_leaders,idx_followers)
    global Npart L linewidth slow_interval
    skip = 2;
    R = L/2*pi;
    
    % Slow interval
    
    theta_slow = 2*pi*slow_interval/L;
    
    theta = linspace(0,2*pi,10000);
    x = R*cos(theta);
    y = R*sin(theta);
    
%     if (length(idx_leaders)>0 && length(idx_leaders)<490)
%         idx_leaders
%         setdiff(1:500,idx_leaders)
%         cccc
%     end
    
    theta_part = xi*2*pi/L;
    x_part = R*cos(theta_part);
    y_part = R*sin(theta_part);
        
    theta_leaders = xi(idx_leaders(1:skip:end))*2*pi/L;
    x_leaders = R*cos(theta_leaders);
    y_leaders = R*sin(theta_leaders);
    
    theta_followers = xi(idx_followers(1:skip:end))*2*pi/L;
    x_followers = R*cos(theta_followers);
    y_followers = R*sin(theta_followers);
    
    f = figure(1);
%     set(f,'units','normalized','outerposition',[0 0 1 1])
    % Plot track
    plot(x,y,'k', 'linewidth',linewidth)
    hold on
    % Plot particles
    plot(x_followers,y_followers,'r.')
%     plot(x_part(1:2:end),y_part(1:2:end),'r.')
    % Plot leaders
    plot(x_leaders,y_leaders,'b.')
    % Plot start and ending particle
    plot(x_part(1),y_part(1),'gs','MarkerFaceColor','g')
    plot(x_part(end),y_part(end),'ms','MarkerFaceColor','m')
    % Plot slow area
    pl1 = line([1.05*R*cos(theta_slow(1)) 0],[1.05*R*sin(theta_slow(1)) 0],'linewidth',2);
    pl2 = line([1.05*R*cos(theta_slow(2)) 0],[1.05*R*sin(theta_slow(2)) 0],'linewidth',2);
    pl1.Color = 'red';
    pl1.LineStyle = '--';
    pl2.Color = 'red';
    pl2.LineStyle = '--';
    
    grid on
    pbaspect([1 1 1])
    hold off
%     axis([78.4 78.6 -0.7 1.2])
end


%%