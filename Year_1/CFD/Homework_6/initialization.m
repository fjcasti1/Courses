function [u,v,Y]=initialization(M,N,hx,hy)
%% Initialize u
u = zeros(M+1,N+2);
% Impose left and right boundaries, only once, since they do not depend on
% the interior nor they are altered.
yju=linspace(-hy/2,2+hy/2,N+2);
for j=1:N+2
    if (yju(j)>=0.5 && yju(j)<=1) % Inlet 1
        u(1,j)=2;
    elseif (yju(j)>=1 && yju(j)<=1.5) % Inlet 2
        u(M+1,j)=-1;
    end
end
%% Initialize v
v = zeros(M+2,N+1);
% Impose top and bottom boundaries, only once, since they do not depend on
% the interior except the outlet Neumann BC, which will be updated after 
% iterations since for the initial condition it is satisfied.
xiv=linspace(-hx/2,4+hx/2,M+2);
for i=1:M+2
    if (xiv(i)>=0.5 && xiv(i)<=1) % Inlet 3
        v(i,N+1)=-1;
    end
end
%% Initialize Y
Y=zeros(M+2,N+2);
yjY=linspace(-hy/2,2+hy/2,N+2);
% Since Y is initialized to zero, the ghost cells are updated except the
% ones corresponding to the inlets 1 and 2.
% Inlet 1
a=find(yjY>=0.5,1,'first');
b=find(yjY<=1,1,'last');
for j=a:b
    Y(1,j)=2-Y(2,j);
end
% Inlet 2
a=find(yjY>=1,1,'first');
b=find(yjY<=1.5,1,'last');
for j=a:b
    Y(M+2,j)=0.5-Y(M+1,j);
end
end
end