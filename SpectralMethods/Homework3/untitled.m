clear variables; close all
N = 20;  
  [D,x] = cheb(N); y = x;
  dt = 6/N^2;
  [xx,yy] = meshgrid(x(2:N),y(2:N));
  plotgap = round((1/3)/dt); dt = (1/3)/plotgap;
  vv = exp(-40*((xx-.4).^2 + yy.^2));
  uu = vv;
% Compute Cheb differentiation matrix
    D2 = D^2; D2 = D2(2:N,2:N);
    I = eye(N-1);
% Time-stepping by leap frog formula:
  [ay,ax] = meshgrid([.56 .06],[.1 .55]); clf
  
  
  for n = 0:3*plotgap
    t = n*dt;
    if rem(n+.5,plotgap)<1     % plots at multiples of t=1/3
      i = n/plotgap+1;
      subplot('position',[ax(i) ay(i) .36 .36])
      [xxx,yyy] = meshgrid(-1:1/16:1,-1:1/16:1);
      vvv = interp2(xx,yy,uu,xxx,yyy,'cubic');
      mesh(xxx,yyy,vvv), axis([-1 1 -1 1 -0.15 1])
      colormap(1e-6*[1 1 1]); title(['t = ' num2str(t)]), drawnow
    end
  tic    
    vv_long = vv(:);
    L = kron(I,D2)+kron(D2,I);
% Reshape long 1D results onto 2D grid:
%     uu = reshape(wave,N-1,N-1);
% Computation Time
    A = [zeros((N-1)^2) eye((N-1)^2); L zeros((N-1)^2)];
    dw = [vv_long; zeros((N-1)^2,1)];
    wt = expm(A*t)*dw;
    u = wt(1:(N-1)^2,:);
    uu = reshape(u,N-1,N-1);
%     vvnew = 2*vv - vvold + dt^2*(uu); 
%     vvold = vv; vv = vvnew;
  end
%   time_mat(k) = toc;