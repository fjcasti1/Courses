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
 
% laplacian
L = kron(eye(N-2),D2(2:end-1,2:end-1))+kron(D2(2:end-1,2:end-1),eye(N-2));  %Laplacian of the interior
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
T(1,:) = TempBC(t,x);

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
   tic 
   wi = w(2:end-1,2:end-1); wi=wi(:);
   psi(2:end-1,2:end-1) = reshape(up\(lo\(-wi(per))),N-2,N-2);
   time1(i) = toc;
   tic
   psi2(2:end-1,2:end-1) = sylvester(D2(2:end-1,2:end-1),D2p(2:end-1,2:end-1),-w(2:end-1,2:end-1));
   time2(i) = toc;
   diff(i) = norm(psi-psi2);

   % Update Velocity
   u = D*psi;
   v = -psi*Dp;
   
   % BC's for u,v,T. Vorticity is calculated from u,v. Stream-function is
   % obtained from w.
   u(indb) = 0;
   v(indb) = 0;
   T(:,1) = 0;
   T(:,end) = 0;
   T(end,:) = 0;
   T(1,:) = TempBC(t,x);
   
   % Advance time
   t = t+dt;
  
   count = count + 1;
   if count == 200
       
%        clev = [-734.8568e-009   -35.1228e-006    -1.5250e-003   -12.5292e-003   -22.7080e-003...
%    -31.9941e-003   -44.0936e-003   -57.0353e-003   -69.0273e-003   -78.5768e-003   -84.4691e-003...
%    -85.9000e-003   -82.4650e-003   -74.1142e-003   -61.1339e-003   -44.1653e-003   -23.8358e-003...
%     -5.3693e-003   494.0942e-006   610.6630e-006   159.8619e-006    13.2590e-006    15.1974e-009];
       subplot(2,2,1)
       contourf(xx,yy,psi)
       axis([0 1 0 1]), axis square
       colormap(hot)
       title('streamlines','fontsize',16)
        
       subplot(2,2,2)
       contourf(xx,yy,w,30)
       axis([0 1 0 1]), axis square
       title('vorticity','fontsize',16)
      
       subplot(2,2,3)
       contourf(xx,yy,T)
       axis([0 1 0 1]), axis square
       title('Temperature field','fontsize',16)
       colormap(jet)
       colorbar
       caxis([0 1])
       
       speed = sqrt(u.^2+v.^2);
       
       subplot(2,2,4)
       quiver(xx,yy,u./speed,v./speed)
       axis([0 1 0 1]), axis square
       title('normalized velocity','fontsize',16)
       
       drawnow
              
       count = 0;
   end
   t
end
%%
n = 1:i;
figure
plot(n,time1./time2)
hold on
plot(n,mean(time1./time2)*ones(size(N)),'r')
plot(n,std(time1./time2)*ones(size(N))+mean(time1./time2)*ones(size(N)),'r-.')
plot(n,-std(time1./time2)*ones(size(N))+mean(time1./time2)*ones(size(N)),'r-.')
grid on
xlim([0 n(end)])
xlabel('$n$','interpreter','latex','fontsize',14)
ylabel('$time_1/time_2$','interpreter','latex','fontsize',14)
ylabel('$t_{LU}/t_{Sylvester}$','interpreter','latex','fontsize',14)
saveas(gcf,'Latex/FIGURES/P5_comparison','png')

figure
plot(n,diff)
grid on
xlabel('$n$','interpreter','latex','fontsize',14)
ylabel('$\|\psi_{LU}-\psi_{Sylvester}\|$','interpreter','latex','fontsize',14)
ylabel('$\|\psi_{LU}-\psi_{Sylvester}\|_2$','interpreter','latex','fontsize',14)
xlim([0 n(end)])
saveas(gcf,'Latex/FIGURES/P5_diff','png')


function Ty0 = TempBC(t,x)
    Ty0 = 2^9*(tanh(100*t))^4*x.^5.*(x-1).^4;
end