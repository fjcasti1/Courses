%% HOMEWORK 3 - FRANCISCO CASTILLO
clear all; close all;
labelfontsize = 14;
figformat = 'png';
% % % %% Problem 1 - 8.4 ATAP
% % % % 
% % % rho = 0.5*pi+sqrt(0.25*pi^2-1);
% % % orderAcuracy('tan(x)',50,2,rho)
% % % 
% % % rho = (0.5*pi+sqrt(0.25*pi^2+1));
% % % orderAcuracy('tanh(x)',40,2,rho)
% % % 
% % % rho = (3+sqrt(8));
% % % orderAcuracy('log((x+3)/4)/(x-1)',30,2,rho,'log1')
% % % 
% % % k = atan(pi/2);
% % % rho = k+sqrt(k^2-1);
% % % orderAcuracy('tan(tan(x))',460,20,rho)
% % % 
% % % rho = -1;
% % % orderAcuracy('(1+x)*log(1+x)',460,20,rho,'(1+x)log(1+x)wrong')
% % % orderAcuracy('(1+x)*log(1+x)',460,20,rho,'(1+x)log(1+x)right')
% % % 
% % % nn = 2:2:20;
% % % err = 0*nn;
% % % for k = 1:length(nn)
% % %     N = nn(k);
% % %     x = -.98:0.02:1;
% % %     F = 0*x;
% % %     FN = F;
% % %     for j = 1:length(x) 
% % %         f = chebfun('cos(t^2)', [-1 x(j)]);
% % %         fN = chebfun('cos(t^2)', [-1 x(j)],N);
% % %         F(j) = sum(f);
% % %         FN(j) = sum(fN);
% % %         F = [0 F]; FN = [0 FN]; x = [-1 x];
% % %     end
% % %     err(k) = norm(F-FN);
% % % end
% % % figure
% % % semilogy(nn,err,'r*')
% % % grid on
% % % xlabel('$N$','interpreter','latex')
% % % ylabel('$\|f-p_N\|$','interpreter','latex')
% % % set(gca,'fontsize',labelfontsize)
% % % txt='Latex/FIGURES/integral';
% % % saveas(gcf,txt,figformat)
% % % 
% % % function orderAcuracy(func,Nmax,Nstep,rho,namefig)
% % %     labelfontsize = 14;
% % %     figformat = 'png';
% % %     if nargin < 5
% % %         namefig = func;
% % %     end
% % %     f = chebfun(func);
% % %     nn = 0:Nstep:Nmax; ee = 0*nn;
% % %     for j=1:length(nn)
% % %         n = nn(j);
% % %         fn = chebfun(f,n+1);
% % %         ee(j) = norm(f-fn);
% % %     end
% % %     figure
% % %     if strcmp(namefig,'(1+x)log(1+x)right')
% % %         semilogy(nn,2000*nn.^(-3.8),'-b')
% % %     else 
% % %         semilogy(nn,rho.^(-nn),'-b')
% % %     end
% % %     hold on
% % %     semilogy(nn,ee,'r*')
% % %     grid on
% % %     xlabel('$N$','interpreter','latex')
% % %     ylabel('$\|f-p_N\|$','interpreter','latex')
% % %     set(gca,'fontsize',labelfontsize)
% % %     txt=['Latex/FIGURES/' namefig];
% % %     saveas(gcf,txt,figformat)
% % % end

%% Problem 2 - 8.7 ATAP

% % % rho=1+sqrt(2);
% % % orderAcuracy('sqrt(x^2+1)',50,2,rho)
% % % 
% % % function orderAcuracy(func,Nmax,Nstep,rho,namefig)
% % %     labelfontsize = 14;
% % %     figformat = 'png';
% % %     if nargin < 5
% % %         namefig = func;
% % %     end
% % %     f = chebfun(func);
% % %     nn = 0:Nstep:Nmax; ee = 0*nn;
% % %     for j=1:length(nn)
% % %         n = nn(j);
% % %         fn = chebfun(f,n+1);
% % %         ee(j) = norm(f-fn);
% % %     end
% % %     figure
% % %     if strcmp(namefig,'(1+x)log(1+x)right')
% % %         semilogy(nn,2000*nn.^(-3.8),'-b')
% % %     else 
% % %         semilogy(nn,rho.^(-nn),'-b')
% % %     end
% % %     hold on
% % %     semilogy(nn,ee,'r*')
% % %     grid on
% % %     xlabel('$N$','interpreter','latex')
% % %     ylabel('$\|f-p_N\|$','interpreter','latex')
% % %     set(gca,'fontsize',labelfontsize)
% % %     txt=['Latex/FIGURES/' namefig];
% % %     saveas(gcf,txt,figformat)
% % % end


% % % %% Problem 3 - 6.7 Trefethen
% % % Nmax = 50; E = zeros(Nmax,1);
% % % for N = 1:Nmax
% % %     [D,x] = cheb(N);
% % %     v = 1./(1+x.^2); vprime = -2*x.*v.^2;    % analytic in [-1,1]
% % %     E(N) = norm(D*v-vprime,inf);
% % % end
% % % % Define ellipse and potential on it.
% % % a = sqrt(2) ; b = 1;
% % % phif = log(0.5*(a+b));
% % % % Plot results:
% % % figure
% % % semilogy(1:Nmax,E(:),'r*')
% % % hold on
% % % semilogy(1:Nmax,exp(-(phif+log(2))*(1:Nmax)))
% % % axis([0 Nmax 1e-16 1e3]), grid on
% % % set(gca,'xtick',0:10:Nmax,'ytick',(10).^(-15:5:0))
% % % xlabel N, ylabel error
% % % txt='Latex/FIGURES/P6_7';
% % % saveas(gcf,txt,figformat)
% % % %% Problem 4 - 8.4 Trefethen
% % % format long
% % % % Grid and initial data:
% % % Nvector = [24,48];
% % % time = zeros(length(Nvector),1);
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
% % %     % -- Using fft -- % Done in Program 20
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
% % %     time(k) = toc;
% % % end

%% Problem 5 - 8.5 Trefethen

% % % % Grid and initial data:
% % % 
% % % Nvector = [24];
% % % time = zeros(length(Nvector),1);
% % % for k = 1:length(Nvector)
% % %     N=Nvector(k);
% % %     [D,x] = cheb(N);  y = x;
% % %     D2 = D^2;
% % %     D2 = D2(2:end-1,2:end-1);
% % %     L = kron(eye(N-1),D2)+kron(D2,eye(N-1));
% % %     A = [zeros(size(L)) eye(size(L)); L zeros(size(L))];
% % %     dt = 6/N^2;
% % %     [xx,yy] = meshgrid(x(2:N),y(2:N));
% % %     x = xx(:); y = yy(:);
% % % 
% % %     plotgap = round((1/3)/dt); dt = (1/3)/plotgap;
% % %     u0 = exp(-40*((x-.4).^2 + y.^2));
% % %     u = u0;
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
% % %     % ------------------------------ %    
% % %     % -- Using matrix exponential -- %
% % %     % ------------------------------ %    
% % %         v = expm(A*t)*[u0;zeros(size(u0))];
% % %         u = v(1:length(u0));
% % %     end
% % %     time(k) = toc;
% % % end

%% Problem 7
close all
f = chebfun('3/(5-4*cos(x))',[-pi,pi]);
plot(f)
grid on
N = 2:2:100;
for k = 1:length(N)
    fcheb = chebfun('3/(5-4*cos(x))',[-pi,pi],N(k));
    ffour = chebfun('3/(5-4*cos(x))',[-pi,pi],N(k),"trig");
    errcheb(k) = norm(f-fcheb,inf);
    errfour(k) = norm(f-ffour,inf);
    % b
    [~,x] = cheb(N(k)); x = pi*x;
    hcheb(k) = max(abs(x(2:end)-x(1:end-1)));
    hfour(k) = 2*pi/(N(k));
end
% a
figure
semilogy(N,errcheb,'b',N,errfour,'r')
hold on
semilogy(N,errcheb,'b*',N,errfour,'r*')
grid on
xlabel('$N$','interpreter','latex')
ylabel('$Error$','interpreter','latex')
set(gca,'fontsize',labelfontsize)
legend('Chebishev', 'Fourier')
txt='Latex/FIGURES/P7_a';
saveas(gcf,txt,figformat)
% b
figure
semilogy(hcheb.^(-1),errcheb,'b',hfour.^(-1),errfour,'r')
hold on
semilogy(hcheb.^(-1),errcheb,'b*',hfour.^(-1),errfour,'r*')
grid on
xlabel('$1/h$','interpreter','latex')
ylabel('$Error$','interpreter','latex')
set(gca,'fontsize',labelfontsize)
legend('Chebishev', 'Fourier')
txt='Latex/FIGURES/P7_b';
saveas(gcf,txt,figformat)
% c
figure
plot(N,hcheb./hfour,'r*')
grid on
axis([0 N(end) 0 pi])
xlabel('$N$','interpreter','latex')
ylabel('$h_{Cheb}/h_{Fourier}$','interpreter','latex')
set(gca,'fontsize',labelfontsize)
legend('Chebishev', 'Fourier')
txt='Latex/FIGURES/P7_c';
saveas(gcf,txt,figformat)
