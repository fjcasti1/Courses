%% Simple code to solve the diffusion equation
% u_t = u_xx -1<x<1, with Neumann boundary conditions

N = 50;
[D,x] = cheb(N);
x = -x;
D = -D;
D2 = D^2;
%%


% initial condition
u0 = cos(pi*x)+0.2*cos(5*pi*x)+1;
% boundary conditions
g1 = @(t) 0*t+1;
g2 = @(t) 0*t+1;

B = D2(2:end-1,[1 end]);
D2 = D2(2:end-1,2:end-1);

A = [D(1,1) D(1,end); D(end,1) D(end,end)];
C = [D(1,2:end-1); D(end,2:end-1)];

t = 0:0.001:2;
%%
[T,U] = ode45(@(t,u) D2*u+(B/A)*([g1(t);g2(t)]-C*u), t, u0(2:end-1));
for k = 2:length(t)
%    u = expm(D2*t(k))*u0(2:end-1);
    plot(x(2:end-1),U(k,:)','*-')
    ylim([-1 2])
    shg
    drawnow
end

    
