%% 2D wave equation Chebyshev+Leap-frog
%% ZERO BC'S
[D,x] = cheb(64);
x = x(2:end-1);
D2 = D^2;
D2 = D2(2:end-1,2:end-1);
 
% 2D grid
[X,Y] = meshgrid(x,x);
 
 
 surf(X,Y,X+Y)
 plot3(X
 cccc
 
% Initial condition
u0 = exp(-20*(X.^2+Y.^2));
u = u0;
 
h = 1-x(2);
dt = h/3.2;
 
% Eigenvalue computations
%L = kron(eye(size(D2)),D2)+kron(D2,eye(size(D2)));
%ee = eig(L);
%plot(real(dt^2*ee),imag(dt^2*ee),'*');
%pause
 
for k = 1:1000
    for j = 1:10
        u2 = 2*u- u0 + dt^2*(D2*u+u*D2');
        u0 = u;
        u = u2;
    end
    surf(X,Y,u)
    zlim([-1 1])
    drawnow
    shg
end