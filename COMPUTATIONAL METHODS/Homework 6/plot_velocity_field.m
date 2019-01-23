% This code plot the velocity filed of a flow in a cavity
% The data is recovered from the file velocitydata.mat
%

load velocitydata.mat

N = size(xx,1);

     % Use a subset of gridpoints
       ip = 1:4:N+1;
       up = u(ip,ip); vp = v(ip,ip);
       xp = xx(ip,ip); yp = yy(ip,ip);
       
       % Velocity
       subplot(1,2,1)
       quiver(xp,yp,up,vp,10)
       axis([0 1 0 1]), axis square
       title('velocity field','fontsize',16)
       
       
       % Stream function
       subplot(1,2,2) 
       mp = max(max(str));
       mm = min(min(str));
       contour(xx,yy,str,-logspace(-10,log10(-mm),50),'r'),  hold on
       contour(xx,yy,str, logspace(-10,log10(mp),20),'k'),   hold off
       axis square, title('streamlines (logscale)','fontsize',16)
       
       