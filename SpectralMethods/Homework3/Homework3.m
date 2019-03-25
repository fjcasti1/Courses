%% HOMEWORK 3 - FRANCISCO CASTILLO
clear all; close all;

% % % %% Problem 5 - 8.4 Trefethen
% % % 
% % % % Grid and initial data:
% % % 
% % % Nvector = [24,48];
% % % time2 = zeros(length(Nvector),1);
% % % for k = 1:length(Nvector)
% % %     N=Nvector(k);
% % %     [D,x] = cheb(N);  y = x;
% % %     D2 = D^2;
% % %     D2 = D2(2:end-1,2:end-1);
% % %     L = kron(eye(N-1),D2)+kron(D2,eye(N-1));
% % %     dt = 6/N^2;
% % %     [xx,yy] = meshgrid(x(2:N),y(2:N));
% % %     x = xx(:); y = yy(:);
% % % 
% % %     plotgap = round((1/3)/dt); dt = (1/3)/plotgap;
% % %     u = exp(-40*((x-.4).^2 + y.^2));
% % %     uold = u; 
% % % 
% % %     % Time-stepping by leap frog formula:
% % %     [ay,ax] = meshgrid([.56 .06],[.1 .55]); clf
% % %     tic
% % %     for n = 0:3*plotgap
% % %         t = n*dt;
% % %         if rem(n+.5,plotgap)<1     % plots at multiples of t=1/3
% % %           uu = reshape(u,N-1,N-1);
% % %           i = n/plotgap+1;
% % %           subplot('position',[ax(i) ay(i) .36 .36])
% % %           [xxx,yyy] = meshgrid(-1:1/16:1,-1:1/16:1);
% % %           uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
% % %           mesh(xxx,yyy,uuu), axis([-1 1 -1 1 -0.15 1])
% % %           colormap(1e-6*[1 1 1]); title(['t = ' num2str(t)]), drawnow
% % %         end
% % %     % -------------------- %    
% % %     % -- Using matrices -- %
% % %     % -------------------- %    
% % %         unew = 2*u - uold + dt^2*L*u; 
% % %         uold = u; u = unew;
% % %     % --------------- %
% % %     % -- Using fft -- %
% % %     % --------------- %
% % %     %     uxx = zeros(N+1,N+1); uyy = zeros(N+1,N+1);
% % %     %     ii = 2:N;
% % %     %     for i = 2:N                % 2nd derivs wrt x in each row
% % %     %       v = vv(i,:); V = [v fliplr(v(ii))];
% % %     %       U = real(fft(V));
% % %     %       W1 = real(ifft(1i*[0:N-1 0 1-N:-1].*U)); % diff wrt theta
% % %     %       W2 = real(ifft(-[0:N 1-N:-1].^2.*U));    % diff^2 wrt theta
% % %     %       uxx(i,ii) = W2(ii)./(1-x(ii).^2) - x(ii).* ... 
% % %     %                      W1(ii)./(1-x(ii).^2).^(3/2);
% % %     %     end
% % %     %     for j = 2:N                % 2nd derivs wrt y in each column
% % %     %       v = vv(:,j); V = [v; flipud(v(ii))];
% % %     %       U = real(fft(V));
% % %     %       W1 = real(ifft(1i*[0:N-1 0 1-N:-1]'.*U));% diff wrt theta   
% % %     %       W2 = real(ifft(-[0:N 1-N:-1]'.^2.*U));   % diff^2 wrt theta
% % %     %       uyy(ii,j) = W2(ii)./(1-y(ii).^2) - y(ii).* ...
% % %     %                      W1(ii)./(1-y(ii).^2).^(3/2);
% % %     %     end
% % %     end
% % %     time2(k) = toc;
% % % end
%% Problem 5 - 8.5 Trefethen

% Grid and initial data:

Nvector = [20,24,48];
time2 = zeros(length(Nvector),1);
for k = 1:1%length(Nvector)
    N=Nvector(k);
    [D,x] = cheb(N);  y = x;
    D2 = D^2;
    D2 = D2(2:end-1,2:end-1);
    L = kron(eye(N-1),D2)+kron(D2,eye(N-1));
    A = [zeros(size(L)) eye(size(L)); L zeros(size(L))];
    dt = 6/N^2;
    [xx,yy] = meshgrid(x(2:N),y(2:N));
    x = xx(:); y = yy(:);

    plotgap = round((1/3)/dt); dt = (1/3)/plotgap;
    u0 = exp(-40*((x-.4).^2 + y.^2));
    u = u0;
    % Time-stepping by leap frog formula:
    [ay,ax] = meshgrid([.56 .06],[.1 .55]); clf
    tic
    for n = 0:3*plotgap
        t = n*dt;
        if rem(n+.5,plotgap)<1     % plots at multiples of t=1/3
          uu = reshape(u,N-1,N-1);
          i = n/plotgap+1;
          subplot('position',[ax(i) ay(i) .36 .36])
          [xxx,yyy] = meshgrid(-1:1/16:1,-1:1/16:1);
          uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
          mesh(xxx,yyy,uuu), axis([-1 1 -1 1 -0.15 1])
          colormap(1e-6*[1 1 1]); title(['t = ' num2str(t)]), drawnow
        end
    % ------------------------------ %    
    % -- Using matrix exponential -- %
    % ------------------------------ %    
        v = expm(A*t)*[u0;zeros(size(u0))];
        u = v(1:length(u0));
    end
    time2(k) = toc;
end