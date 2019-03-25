function [u,x]=PDE_solve(N,tmax,eps,movie)
x0 = 0;
xf = 2*pi;
h = (xf-x0)/N;
x = x0 + h*(1:N)';
u0 = sin(x/2).^4;
 
tic;
[t,u] = ode113(@(t,u) rhs(u,N,eps), [0 tmax], u0);
toc
if movie
    figure
    for k = 1:2:length(t)
        plot(x,u(k,:));
        axis([x0 xf -.2 1.2]);
        drawnow
    end
end

end
 
function y = rhs(u,N,eps)
    vhat = fft(u);
    what = 1i*[0:N/2-1 0 -N/2+1:-1]' .* vhat;
    what2 = -([0:N/2-1 N/2 -N/2+1:-1].^2)' .* vhat;
    
    ux = real(ifft(what));
    uxx = real(ifft(what2));
    
    y = 4*pi^2*eps*uxx-2*pi*u.*ux;
end