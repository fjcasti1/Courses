function [pu,pv,pY]=probeValues(u,v,Y,M,N,hx,hy)
% Define the points of the different meshes
xu=linspace(0,4,M+1);
yu=linspace(-hy/2,2+hy/2,N+2);
xv=linspace(-hx/2,4+hx/2,M+2);
yv=linspace(0,2,N+1);
xY=linspace(-5*hx/2,4+5*hx/2,M+6);
yY=linspace(-5*hx/2,2+5*hx/2,N+6);
% Probe for u
pu = ( u(xu==1,find(yu<=0.5,1,'last')) + u(xu==1,find(yu>=0.5,1)) )/2;
% Probe for v
pv = ( v(find(xv<=1,1,'last'),yv==1.5) + v(find(xv>=1,1),yv==1.5) )/2;
% Probe for Y
pY = ( Y(find(xY<=2,1,'last'),find(yY<=0.5,1,'last'))...
    +Y(find(xY<=2,1,'last'),find(yY>0.5,1))...
    +Y(find(xY>2,1),find(yY<=0.5,1,'last'))...
    +Y(find(xY>2,1),find(yY>0.5,1)) )/4;
end