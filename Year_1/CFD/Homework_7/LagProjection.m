function [u,v]=LagProjection(u,v,phi,dt,hx,hy)
M=size(phi,1)-2;
N=size(phi,2)-2;
for i=1:M+1
    for j=1:N+2
        phix(i,j)=(phi(i+1,j)-phi(i,j))/hx;
    end
end
for i=1:M+2
    for j=1:N+1
        phiy(i,j)=(phi(i,j+1)-phi(i,j))/hy;
    end
end
u=u-dt*phix;
v=v-dt*phiy;
% Apply boundary conditions to ghost cells only
u(:,1)=-u(:,2);
u(:,N+2)=-u(:,N+1);
v(1,:)=-v(2,:);
v(M+2,:)=-v(M+1,:);
end