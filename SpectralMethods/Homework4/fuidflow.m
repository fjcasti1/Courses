%% Lid Driven Cavity - Rodrigo Platte
clear all; close all; clc;

% Parameters
N = 40;  %40
Re = 2000;
dt = 1e-4;
 
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
u = 0*xx;
v = u;
str = u; 
% Everything initialized to zero
 
% boundary condition
u(end,:) = 1;%16*(x.*(1-x)).^2; % END ??? check xx yy and see differently
ub = u;
 
count = 0;
t = 0;
% main loop
while t<200
    
   % vorticity 
   w = v*Dp-D*u;  % w = dvdx-dudy
    
   % advance vorticity
   w = w + dt*(-u.*(w*Dp)-v.*(D*w)+(1/Re)*(w*D2p+D2*w));
    
   % compute stream function
   wi = w(2:end-1,2:end-1); wi=wi(:);
   str(2:end-1,2:end-1) = reshape(up\(lo\(-wi(per))),N-2,N-2);
    
   % update velocity
   u = D*str;
   v = -str*Dp;
    
   % BC's for u and v
   u(indb) = ub(indb);
   v(indb) = 0;
    
   t = t+dt;
   count = count + 1;
   if count == 100
        
       subplot(2,2,1)
       quiver(xx,yy,u,v,10)
       axis([0 1 0 1]), axis square
       title('velocity field','fontsize',16)
        
       subplot(2,2,2)
       speed = sqrt(u.^2+v.^2);
       quiver(xx,yy,u./speed,v./speed)
       axis([0 1 0 1]), axis square
       title('normalized velocity','fontsize',16)
        
       subplot(2,2,3) 
       clev = [-734.8568e-009   -35.1228e-006    -1.5250e-003   -12.5292e-003   -22.7080e-003...
   -31.9941e-003   -44.0936e-003   -57.0353e-003   -69.0273e-003   -78.5768e-003   -84.4691e-003...
   -85.9000e-003   -82.4650e-003   -74.1142e-003   -61.1339e-003   -44.1653e-003   -23.8358e-003...
    -5.3693e-003   494.0942e-006   610.6630e-006   159.8619e-006    13.2590e-006    15.1974e-009];
       contourf(xx,yy,str,clev)
       axis([0 1 0 1]), axis square
       colormap(hot)
       title('streamlines','fontsize',16)
        
       subplot(2,2,4)
       w = v*Dp-D*u;
       contourf(xx,yy,w,30)
       axis([0 1 0 1]), axis square
       title('vorticity','fontsize',16), drawnow
        
       count = 0;
        
   end
 
end