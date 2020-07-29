%% Homework 4, Problem 5 -  Francisco Castillo
clear all; close all; clc;

% Parameters
N = 40;  %40
Pr = 0.71;
Ra = 2e5;
dt = 2e-6;
 
% Grid and diff matrices
[D,xch] = cheb(N-1);
x = (-xch+1)/2; D = -2*D;   % To trasnlate the domain to [0,1]
[xx,yy] = meshgrid(x);

Dp = D';
D2 = D^2; D2p = D2';
indb = find(xx==0|xx==1|yy==0|yy==1); % For what?? Impose BCs in velocities
 

% Laplacian of the interior
L = kron(eye(N-2),D2(2:end-1,2:end-1))+kron(D2(2:end-1,2:end-1),eye(N-2));  
%Linv = inv(L);
[lo,up,per] = lu(L,'vector');  % LU factorization
 
% Initial velocity & pre-allocate memory
T = 0*xx;
psi = T;
psi2 = T;
w = T;
u = T;
v = T;
% Everything initialized to zero
count = 0;
t = 0;

% boundary condition
Topt = 'natural convection';
if (strcmp(Topt,'natural convection'))
    T(1,:) = TempBC(t,x);
elseif (strcmp(Topt,'blocked convection'))
	T(end,:) = TempBC(t,x);
end

% main loop
i=0;
while t<.1%200
    i=i+1;
    % vorticity 
    w = v*Dp-D*u;  % w = dvdx-dudy

    % Advance Vorticity
    w =  w + dt*(-v.*(D*w)-u.*(w*Dp)+Pr*(w*D2p+D2*w)+Ra*Pr*T*Dp);

    % Advance Temperature
    T = T + dt*(-v.*(D*T)-u.*(T*Dp)+T*D2p+D2*T);

    % compute stream function
    %    tic 
    %    wi = w(2:end-1,2:end-1); wi=wi(:);
    %    psi(2:end-1,2:end-1) = reshape(up\(lo\(-wi(per))),N-2,N-2);
    %    time1(i) = toc;
    %    tic
    psi(2:end-1,2:end-1) = sylvester(D2(2:end-1,2:end-1),D2p(2:end-1,2:end-1),-w(2:end-1,2:end-1));
    %    time2(i) = toc;
    %    diff(i) = norm(psi-psi2);

    % Update Velocity
    u = D*psi;
    v = -psi*Dp;

    % BC's for u,v,T. Vorticity is calculated from u,v. Stream-function is
    % obtained from w.
    u(indb) = 0;
    v(indb) = 0;
    T(:,1) = 0;
    T(:,end) = 0;
    if (strcmp(Topt,'natural convection'))
        T(1,:) = TempBC(t,x); 
        T(end,:) = 0;
    elseif (strcmp(Topt,'blocked convection'))
        T(1,:) = 0;
        T(end,:) = TempBC(t,x);
    end

    % Advance time
    t = t+dt;

   count = count + 1;
   if count == 200
       
       subplot(2,2,1)
       contourf(xx,yy,psi)
       axis([0 1 0 1]), axis square
       colormap(hot)
       title('Streamlines','fontsize',16)
        
       subplot(2,2,2)
       contourf(xx,yy,w,30)
       axis([0 1 0 1]), axis square
       title('Vorticity','fontsize',16)
      
       subplot(2,2,3)
       contourf(xx,yy,T)
       axis([0 1 0 1]), axis square
       title('Temperature','fontsize',16)
       colormap(jet)
       colorbar
       caxis([0 1])
       
       speed = sqrt(u.^2+v.^2);
       
       subplot(2,2,4)
       quiver(xx,yy,u./speed,v./speed)
       axis([0 1 0 1]), axis square
       title('Velocity','fontsize',16)
       
       drawnow
              
       count = 0;
   end

end
if (strcmp(Topt,'natural convection'))
    saveas(gcf,'Latex/FIGURES/P5','png')
elseif (strcmp(Topt,'blocked convection'))
    saveas(gcf,'Latex/FIGURES/P5_blocked','png')
end

%%
% n = 1:i;
% figure
% plot(n,time1./time2)
% hold on
% plot(n,mean(time1./time2)*ones(size(n)),'r')
% plot(n,(mean(time1./time2)+std(time1./time2))*ones(size(n)),'r-.')
% plot(n,(mean(time1./time2)-std(time1./time2))*ones(size(n)),'r-.')
% grid on
% xlim([0 n(end)])
% xlabel('$n$','interpreter','latex','fontsize',14)
% ylabel('$time_1/time_2$','interpreter','latex','fontsize',14)
% ylabel('$t_{LU}/t_{Sylvester}$','interpreter','latex','fontsize',14)
% saveas(gcf,'Latex/FIGURES/P5_comparison','png')
% 
% figure
% plot(n,diff)
% grid on
% xlabel('$n$','interpreter','latex','fontsize',14)
% ylabel('$\|\psi_{LU}-\psi_{Sylvester}\|$','interpreter','latex','fontsize',14)
% ylabel('$\|\psi_{LU}-\psi_{Sylvester}\|_2$','interpreter','latex','fontsize',14)
% xlim([0 n(end)])
% saveas(gcf,'Latex/FIGURES/P5_diff','png')


function Ty0 = TempBC(t,x)
    Ty0 = 2^9*(tanh(100*t))^4*x.^5.*(x-1).^4;
end