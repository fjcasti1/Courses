function [u,v,Y]=initialization(M,N)
t=0;
%% Initialize u
u = zeros(M+1,N+2);
% Boundary Conditions
u=applyBCs(u,M,N,t,'u');

%% Initialize v
v = zeros(M+2,N+1);
% Boundary Conditions
v=applyBCs(v,M,N,t,'v');

%% Initialize Y
Y = zeros(M+6,N+6);
% Boundary Conditions
Y=applyBCs(Y,M,N,t,'Y');
end