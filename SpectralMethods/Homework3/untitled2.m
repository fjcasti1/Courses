%% Problem 8.4
%% 
% 2nd-order wave eq. in 2D via Cheb matrices
clear variables;
% Grid and initial data:
  gridsize = [24,48];
for k = 1:length(gridsize)
  N = gridsize(k);  
  [D,x] = cheb(N); y = x;
  dt = 6/N^2;
  [xx,yy] = meshgrid(x(2:N),y(2:N));
  plotgap = round((1/3)/dt); dt = (1/3)/plotgap;
  vv = exp(-40*((xx-.4).^2 + yy.^2));
  vvold = vv; 
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
      vvv = interp2(xx,yy,vv,xxx,yyy,'cubic');
      mesh(xxx,yyy,vvv), axis([-1 1 -1 1 -0.15 1])
      colormap(1e-6*[1 1 1]); title(['t = ' num2str(t)]), drawnow
    end
  tic    
    vv_long = vv(:);
    uxx = kron(I,D2)*vv_long;
    uyy = kron(D2,I)*vv_long;    
    wave = uxx+uyy;
    
% Reshape long 1D results onto 2D grid:
    uu = reshape(wave,N-1,N-1);
% Computation Time
    vvnew = 2*vv - vvold + dt^2*(uu); 
    vvold = vv; vv = vvnew;
  end
  time_mat(k) = toc;
end

%% 
% p20.m - 2nd-order wave eq. in 2D via FFT (compare p19.m)
  gridsize = [24,48];
for k = 1:length(gridsize)
  N = gridsize(k);  
% Grid and initial data:
    x = cos(pi*(0:N)/N); y = x';
  dt = 6/N^2;
  [xx,yy] = meshgrid(x,y);
  plotgap = round((1/3)/dt); dt = (1/3)/plotgap;
  vv = exp(-40*((xx-.4).^2 + yy.^2));
  vvold = vv; 

% Time-stepping by leap frog formula:
  [ay,ax] = meshgrid([.56 .06],[.1 .55]); clf
  for n = 0:3*plotgap
    t = n*dt;
    if rem(n+.5,plotgap)<1     % plots at multiples of t=1/3
      i = n/plotgap+1;
      subplot('position',[ax(i) ay(i) .36 .36])
      [xxx,yyy] = meshgrid(-1:1/16:1,-1:1/16:1);
      vvv = interp2(xx,yy,vv,xxx,yyy,'cubic');
      mesh(xxx,yyy,vvv), axis([-1 1 -1 1 -0.15 1])
      colormap(1e-6*[1 1 1]); title(['t = ' num2str(t)]), drawnow
    end
  tic
    uxx = zeros(N+1,N+1); uyy = zeros(N+1,N+1);
    ii = 2:N;
    for i = 2:N                % 2nd derivs wrt x in each row
      v = vv(i,:); V = [v fliplr(v(ii))];
      U = real(fft(V));
      W1 = real(ifft(1i*[0:N-1 0 1-N:-1].*U)); % diff wrt theta
      W2 = real(ifft(-[0:N 1-N:-1].^2.*U));    % diff^2 wrt theta
      uxx(i,ii) = W2(ii)./(1-x(ii).^2) - x(ii).* ... 
                     W1(ii)./(1-x(ii).^2).^(3/2);
    end
    for j = 2:N                % 2nd derivs wrt y in each column
      v = vv(:,j); V = [v; flipud(v(ii))];
      U = real(fft(V));
      W1 = real(ifft(1i*[0:N-1 0 1-N:-1]'.*U));% diff wrt theta   
      W2 = real(ifft(-[0:N 1-N:-1]'.^2.*U));   % diff^2 wrt theta
      uyy(ii,j) = W2(ii)./(1-y(ii).^2) - y(ii).* ...
                     W1(ii)./(1-y(ii).^2).^(3/2);
    end
    % Computation Time
    vvnew = 2*vv - vvold + dt^2*(uxx+uyy); 
    vvold = vv; vv = vvnew;
  end
  time_fft(k) = toc;
end
 
% % % %% Problem 8.6
% % % % chebfft2
% % % clear variables; close all
% % % 
% % % N = 20;
% % % [D,x] = cheb(N);
% % % % example funcitons
% % % 
% % % % f = 50*cos(2*x)+27*sin(10*x);
% % % f = exp(x).*sin(5*x);
% % % 
% % % % compute second derivative using fft
% % % ddf_fft = chebfft(chebfft(f));
% % % % compute derivative using matrices
% % % D2 = D^2;
% % % ddf_mat = (D2)*f;
% % % 
% % % % difference
% % % diff = ddf_fft - ddf_mat;
% % % maxdiff = norm(ddf_fft - ddf_mat,inf);
% % % [largest,I] = max(diff);
% % % 
% % % % Plot the solutions overlapping
% % % figure
% % % subplot(1,2,1)
% % % plot(x,ddf_fft,'r-*',x,ddf_mat,'b--o','LineWidth',2)
% % % grid on
% % % legend('FFT','Matrices')
% % % xlabel('x')
% % % ylabel('Second Derivative')
% % % subplot(1,2,2)
% % % plot(x,diff,'-.','LineWidth',3)
% % % hold on
% % % plot(x(I),largest,'kp','MarkerSize',12,'MarkerFaceColor','k')
% % % legend('Difference in values','largest difference')
% % % xlabel('x')
% % % ylabel('error')
% % % grid on

%% Problem 8.5
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
    whos dw wt u uu
    ccc
%     vvnew = 2*vv - vvold + dt^2*(uu); 
%     vvold = vv; vv = vvnew;
  end
%   time_mat(k) = toc;
cccc
%% Problem 7
N = 20;
h = 2*pi/N;
x = -pi+h*(1:N)';
xx = linspace(-pi,pi,50);
% theta = (0:N)*pi/N;
% x = cos(theta);
f = 3./(5-4*cos(x));

% Using ffts
coeff = fft(f);
coeff = fftshift(coeff/N);
% for j = 1:length(xx)
%    f_fft(j) = coeff.'*exp(1i*(-N/2+1:N/2)*xx(j)).'; 
% end
f_fft = sum(coeff.'.*exp(xx*1i.*(-N:N)'));

% f_fft = coeff'*exp(1i*[0:N/2-1 0 -N/2+1:-1]');
% f_fft = [0:N/2-1 0 -N/2+1:-1]'.*coeff;

% Polynomial Approximation
% f_poly = chebfun


figure
% plot(x,f)
% hold on
plot(x,f_fft)
% plot(x,c)
% legend
% hold off