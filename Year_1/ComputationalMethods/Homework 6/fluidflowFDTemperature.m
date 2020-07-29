% Fluid flow in a cavity (Navier-Stokes equations)
% Stream-function and vorticity formulation
%
% Rodrigo Platte, Arizona State University, April 2013.

function fluidflowFDTemperature(xx,yy,u,v)
% Parameters
path='Latex/FIGURES/';
N = 150;
alpha = 0.0005;
dt = min(0.2/(alpha*N^2),1e-3);

h = 1/N; % N+1 points in each direction

% Initial condition & pre-allocate memory (start with zeros at t = 0)
T = 0*xx;

% boundary condition
ii = 2:N; jj = 2:N;  % index for interior nodes
Ti = zeros(N-1,N-1); % interior values of vorticity

count = 0;
time = 0;
outputTime=[10 50 150 300];
endtime=outputTime(end);
n=1;
% main loop
while time<endtime
	if ( time < outputTime(n) && time+dt >= outputTime(n) ) %#ok<ALIGN>
        dt=outputTime(n)-time;
        n=n+1;
    else
        dt = min(0.2/(alpha*N^2),1e-3);
    end
   % advance Temperature with forward Euler
   Ti = Ti + dt*Trhs(u,v,T,alpha);
   T(ii,jj)=Ti;
   
   % update Temperature at boundary
   T = BCT(xx,yy,time,T);
   
   time = time+dt;
   count = count + 1;
   
   % Plot results every now and then
   if ismember(time,outputTime);
       plotresults(time,T);    
   end

end
    % Right-hand side of vorticity equation
    function rhs = Trhs(u,v,T,alpha)
        Tx = (T(ii,jj+1)-T(ii,jj-1))/(2*h);
        Ty = (T(ii+1,jj)-T(ii-1,jj))/(2*h);
        LT = (T(ii,jj+1)+T(ii,jj-1)+T(ii+1,jj)+T(ii-1,jj)-4*T(ii,jj))/h^2;
        rhs = -u(ii,jj).*Tx - v(ii,jj).*Ty + alpha*LT;          
    end

    % Enforce boundary condition (temperature)
    function T = BCT(xx,yy,t,T)
        T(1,:)=xx(1,:).*yy(1,:).*tanh(10*t).^2;
        T(N+1,:)=xx(N+1,:).*yy(N+1,:).*tanh(10*t).^2;
        T(:,1)=xx(:,1).*yy(:,1).*tanh(10*t).^2;
        T(:,N+1)=xx(:,N+1).*yy(:,N+1).*tanh(10*t).^2;
    end

    % Plot results
    function plotresults(time,T)
        figure();set(gcf,'Visible', 'off');
        pcolor(xx,yy,T)
        shading interp
        H = colorbar; set(H,'fontsize',16)
        caxis([0 1])
        colorbar
        axis square
        drawnow      
        xlabel('$x$','Interpreter','latex')
        ylabel('$y$','Interpreter','latex')
        set(get(gca,'ylabel'),'rotation',0)
        txt=[path,'T_' num2str(time)];
        saveas(gcf,txt,'epsc')
    end

end
    