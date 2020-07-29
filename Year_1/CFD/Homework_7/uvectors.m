function [a,b,c,d] = uvectors(u,M,N,dt,hx,hy,alpha,step)
dx=alpha*dt/hx^2;
dy=alpha*dt/hy^2;
d1=dx/2;
d2=dy/2;
if step==1
    a = -d1*ones((M-1)*N,1);
    b = (1+2*d1)*ones((M-1)*N,1);
    c = -d1*ones((M-1)*N,1);
    d = zeros((M-1)*N,1);
    for j=2:N+1
        for i=2:M
            if i==2 % Left boundary BCs
                a((j-2)*(M-1)+1)=0;
                d((j-2)*(M-1)+1)=d2*u(2,j-1)+(1-2*d2)*u(2,j)+d2*u(2,j+1)...
                                    +d1*u(1,j); %u(1,j) fixed at initialization
            elseif i==M % Right boundary BCs
                c((j-2)*(M-1)+M-1)=0;
                d((j-2)*(M-1)+M-1)=d2*u(i,j-1)+(1-2*d2)*u(i,j)+d2*u(i,j+1)...
                                    +d1*u(M+1,j); %u(M+1,j) fixed at initialization
            else % The top and bottom boundaries are imposed by updating the ghost
                    % cells accordingly with them.
                d((j-2)*(M-1)+i-1)=d2*u(i,j-1)+(1-2*d2)*u(i,j)+d2*u(i,j+1);
            end
        end
    end
elseif step==2
    a = -d2*ones((M-1)*N,1);
    b = (1+2*d2)*ones((M-1)*N,1);
    c = -d2*ones((M-1)*N,1);
    d = zeros((M-1)*N,1);
    for i=2:M
        for j=2:N+1
            if j==2 % Bottom boundary BCs
                a((i-2)*N+1)=0;
                b((i-2)*N+1)=1+3*d2;
                d((i-2)*N+1)=d1*u(i-1,2)+(1-2*d1)*u(i,2)+d1*u(i+1,2);
            elseif j==N+1 % Top boundary BCs
                b((i-1)*N)=1+3*d2;
                c((i-1)*N)=0;
                d((i-1)*N)=d1*u(i-1,N+1)+(1-2*d1)*u(i,N+1)+d1*u(i+1,N+1);
            else % The left and right boundaries are imposed by initialization,
                    % they do not change
                d((i-2)*N+j-1)=d1*u(i-1,j)+(1-2*d1)*u(i,j)+d1*u(i+1,j);
            end
        end
    end
end
end