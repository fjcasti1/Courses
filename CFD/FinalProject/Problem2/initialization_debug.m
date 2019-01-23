function [u,v,Y]=initialization_debug(M,N,hx,hy)
x=linspace(-hx/2,4+hx/2,M+2);
y=linspace(-hy/2,2+hy/2,N+2);
%% Initialize u
u = zeros(M+1,N+2);
for i=2:M
    for j=1:N+2
        u(i,j)=sin(1.2*(x(i)+hx/2)+0.1)+sin(1.4*y(j));
    end
end
% Boundary Conditions
u=applyBCs(u,M,N,hx,hy,'u');

%% Initialize v
v = zeros(M+2,N+1);
for i=1:M+2
    for j=2:N
        v(i,j)=cos(1.4*(y(j)+hy/2))+cos(1.3*x(i));
    end
end
% Boundary Conditions
v=applyBCs(v,M,N,hx,hy,'v');

%% Initialize Y
x=linspace(-5*hx/2,4+5*hx/2,M+6);
y=linspace(-5*hy/2,2+5*hy/2,N+6);
Y = zeros(M+6,N+6);
for i=1:M+6
    for j=1:N+6
        Y(i,j)=sin(0.7*x(i))+cos(1.2*y(j));
    end
end
% Boundary Conditions
Y=applyBCs(Y,M,N,hx,hy,'Y');
end