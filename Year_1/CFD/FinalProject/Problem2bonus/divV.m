function rhs=divV(u,v,M,N,dt)
global hx; global hy;
rhs=zeros(M+2,N+2);
for i=2:M+1
    for j=2:N+1
        rhs(i,j)=(u(i,j)-u(i-1,j))/(hx*dt)+(v(i,j)-v(i,j-1))/(hy*dt);
    end
end
end