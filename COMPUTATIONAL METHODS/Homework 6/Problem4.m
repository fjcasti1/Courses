% Fluid flow in a cavity (Navier-Stokes equations)
% Stream-function and vorticity formulation
%
% Rodrigo Platte, Arizona State University, April 2013.

function Problem4(Ra,Pr)
% Parameters
path='Latex/FIGURES/';
N = 150;
dt = 1e-5;

h = 1/N; % N+1 points in each direction
[xx,yy] = meshgrid(0:h:1); % 2D gridpoints

% Laplacian with zero BCs
o = ones(N-1,1);
D2 = (diag(-2*o) + diag(o(1:N-2),-1)+ diag(o(1:N-2),1))/h^2;
% We need evals and evects to solve Poisson equation
[eV,eval] = eig(D2); lam = diag(eval); 

% Initial condition & pre-allocate memory (start with zeros at t = 0)
u = 0*xx; v = u; w = u; str = u; T=u;

% boundary condition
ii = 2:N; jj = 2:N;  % index for interior nodes
wi = zeros(N-1,N-1); % interior values of vorticity
Ti = zeros(N-1,N-1); % interior values of temperature

count = 0;
time = 0;
outputTime=[1];
endtime=outputTime(end);
n=1;
% main loop
while time<endtime
	if ( time < outputTime(n) && time+dt >= outputTime(n) ) %#ok<ALIGN>
        dt=outputTime(n)-time;
        n=n+1;
    else
        dt = 1e-5;
    end
   % update vorticity at boundary
   w = BCw(u,v,w);
   % update Temperature at boundary
   T = BCT(xx,time,T);
   
   % advance vorticity and temperature with forward Euler
    [rhsT,Tx]=Trhs(u,v,T);
    Ti = Ti + dt*rhsT;
    T(ii,jj)=Ti;
    wi = wi + dt*wrhs(u,v,w,Pr,Ra,Tx);
    w(ii,jj) = wi;
   
   % compute stream function with conjugate gradient
   str(2:end-1,2:end-1) = SolvePoisson(wi);
   
   % update velocity (interior nodes only)
   u(ii,jj) =  (str(ii+1,jj)-str(ii-1,jj))/(2*h);
   v(ii,jj) = -(str(ii,jj+1)-str(ii,jj-1))/(2*h);
   
   time = time+dt;
   count = count + 1;
   
   % Plot results every now and then
   if ismember(time,outputTime);
       plotresults(u,v,str,w,T);    
   end

end

    % Solve Poisson equation (Sylvester equation)
    function S = SolvePoisson(wi)
        ff = -eV'*wi*eV;
        S = ff;
        for j = 1:N-1
            S(:,j) = ff(:,j)./(lam(j)+lam);
        end
        S = eV*S*eV';
    end

    % Right-hand side of vorticity equation
    function rhs = wrhs(u,v,w,Pr,Ra,Tx)
        wx = (w(ii,jj+1)-w(ii,jj-1))/(2*h);
        wy = (w(ii+1,jj)-w(ii-1,jj))/(2*h);
        Lw = (w(ii,jj+1)+w(ii,jj-1)+w(ii+1,jj)+w(ii-1,jj)-4*w(ii,jj))/h^2;
        rhs = -u(ii,jj).*wx - v(ii,jj).*wy + Pr*Lw+Ra*Pr*Tx;          
    end
    % Right-hand side of temperature equation
    function [rhs, Tx] = Trhs(u,v,T)
        Tx = (T(ii,jj+1)-T(ii,jj-1))/(2*h);
        Ty = (T(ii+1,jj)-T(ii-1,jj))/(2*h);
        LT = (T(ii,jj+1)+T(ii,jj-1)+T(ii+1,jj)+T(ii-1,jj)-4*T(ii,jj))/h^2;
        rhs = -u(ii,jj).*Tx - v(ii,jj).*Ty + LT;          
    end

    % Enforce boundary condition (vorticity)
    function w = BCw(u,v,w)
        % three point sided FD formula (second order accurate)
        w(1,:) = -(-3*u(1,:)+4*u(2,:)-u(3,:))/(2*h);      % = -u_y @ y = 0
        w(N+1,:) = (-3*u(N+1,:)+4*u(N,:)-u(N-1,:))/(2*h); % = -u_y @ y = 1
        w(:,1) = (-3*v(:,1)+4*v(:,2)-v(:,3))/(2*h);       % =  v_x @ x = 0
        w(:,N+1) = -(-3*v(:,N+1)+4*v(:,N)-v(:,N-1))/(2*h);% =  v_x @ x = 1
    end
    % Enforce boundary condition (temperature)
    function T = BCT(xx,t,T)
        T(1,:)=2^9*tanh(100*t)^4*xx(1,:).^5.*(xx(1,:)-1).^4;
    end

    % Plot results
    function plotresults(u,v,str,w,T)
        % Use a subset of gridpoints
       ip = 1:3:N+1;
       up = u(ip,ip); vp = v(ip,ip);
       xp = xx(ip,ip); yp = yy(ip,ip);
       
       % Velocity
       subplot(2,2,1)
       quiver(xp,yp,up,vp,10)
       axis([0 1 0 1]), axis square
       title('velocity field','fontsize',16)
       
       % Normalized velocity
       subplot(2,2,2)
       contourf(xx,yy,T)
       axis([0 1 0 1]), axis square
       title('normalized velocity','fontsize',16)
       
       % Stream function
       subplot(2,2,3) 
       mp = max(max(str));
       mm = min(min(str));
       contour(xx,yy,str,-logspace(-10,log10(-mm),30),'r'),  hold on
       contour(xx,yy,str, logspace(-10,log10(mp),20),'k'),   hold off
       axis square, title('streamlines (logscale)','fontsize',16)
       
       % Vorticity
       subplot(2,2,4)
       logw = log10(abs(w));    logw(abs(w)==0) = nan;
       contourf(xx,yy,logw,-3:0.5:2), 
       H = colorbar; set(H,'fontsize',16)
       set(gca,'position',[0.5703 0.1100 0.3347 0.3412])
       axis square, title('vorticity (log(abs(w))','fontsize',16) 
       drawnow
        txt=[path,'solutionProblem4'];
        saveas(gcf,txt,'epsc')
    end

end
    