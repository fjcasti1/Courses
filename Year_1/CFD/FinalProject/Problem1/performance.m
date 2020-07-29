function [R,S]=performance(Y,v,M,N)
global Lx;
global Ly;
global hx;
global hy;
global xv;
% Calculate R
R=sum(sum(Y(4:M+3,4:N+3).*(1-Y(4:M+3,4:N+3))))*hx*hy/(Ly*Lx);
% Calculate S
Lout=1;
% Mesh points for v in x direction
a=find(xv>=1.5,1);
b=find(xv<2.5,1,'last');
% Y values at the bottom boundary
Y_b=(Y(3:M+4,4)+Y(3:M+4,3))/2;
S=(1/Lout)*sum(-v(a:b,1).*Y_b(a:b).*(1-Y_b(a:b)))*hx;
end