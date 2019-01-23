function [u,v,Y]=initialization(M,N,hx,hy)
x=linspace(-hx/2,4+hx/2,M+2);
y=linspace(-hy/2,2+hy/2,N+2);
%% Initialize u
u = zeros(M+1,N+2);
% Boundary Conditions
u=applyBCs(u,M,N,hx,hy,'u');

%% Initialize v
v = zeros(M+2,N+1);
% Boundary Conditions
v=applyBCs(v,M,N,hx,hy,'v');

%% Initialize Y
x=linspace(-5*hx/2,4+5*hx/2,M+6);
y=linspace(-5*hy/2,2+5*hy/2,N+6);
Y = zeros(M+6,N+6);
% Boundary Conditions
Y=applyBCs(Y,M,N,hx,hy,'Y');
end