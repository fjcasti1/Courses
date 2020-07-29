function [Hu,Hv]=HyperbolicBurgers2D(u,v,M,N)
global hx;
global hy;

Hu=zeros(size(u));
Hv=zeros(size(v));
for i=2:M
    for j=2:N+1
        Hu(i,j)=-(u(i+1,j)^2+2*u(i,j)*(u(i+1,j)-u(i-1,j))-u(i-1,j)^2)/(4*hx)...
            -((u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j))-(u(i,j-1)+u(i,j))*(v(i,j-1)+v(i+1,j-1)))/(4*hy);
    end
end
for i=2:M+1
    for j=2:N
        Hv(i,j)=-(v(i,j+1)^2+2*v(i,j)*(v(i,j+1)-v(i,j-1))-v(i,j-1)^2)/(4*hx)...
            -((u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)))/(4*hy);
    end
end
end