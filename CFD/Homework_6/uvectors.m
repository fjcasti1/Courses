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
    % Case j=2 and i=2, includes BCs
    a(1)=0;
    d(1)=d2*u(2,3)+(1-2*d2)*u(2,2)+d1*u(1,2)+d2*u(2,1);
    %Case j=2 and i between [3,M-1], includes BCs
    for i=3:M-1
        d(i-1)=d2*u(i,3)+(1-2*d2)*u(i,2)+d2*u(i,1);
    end
    % Case j=2 and i=M, includes BCs
    c(M-1)=0;
    d(M-1)=d2*u(M,3)+(1-2*d2)*u(M,2)+d1*u(M+1,2)+d2*u(M,1);
    % Case j between [3,N] and i=2, includes BCs
    for j=3:N
        a((j-2)*(M-1)+1)=0;
        d((j-2)*(M-1)+1)=d2*u(2,j+1)+(1-2*d2)*u(2,j)...
            +d2*u(2,j-1)+d1*u(1,j);
    end
    % Interior of the interior, j in [3,N] and i in [3,M-1], no BCs
    for j=3:N
        for i=3:M-1
            d((j-2)*(M-1)+i-1)=d2*u(i,j+1)+(1-2*d2)*u(i,j)+d2*u(i,j-1);
        end
    end
    % Case j between [3,N] and i=M, includes BCs
    for j=3:N
        c((j-2)*(M-1)+M-1)=0;
        d((j-2)*(M-1)+M-1)=d2*u(M,j+1)+(1-2*d2)*u(M,j)+d2*u(M,j-1)...
            +d1*u(M+1,j);
    end
    % Case j=N+1 and i=2, includes BCs
    a((N-1)*(M-1)+1)=0;
    d((N-1)*(M-1)+1)=(1-2*d2)*u(2,N+1)+d2*u(2,N)...
        +d1*u(1,N+1)+d2*u(2,N+2);
    % Case j=N+1 and i [3,M-1], includes BCs
    for i=3:M-1
        d((N-1)*(M-1)+i-1)=(1-2*d2)*u(i,N+1)+d2*u(i,N)...
            +d2*u(i,N+2);
    end
    % Case j=N+1 and i=M, includes BCs
    c((M-1)*N)=0;
    d((M-1)*N)=(1-2*d2)*u(M,N+1)+d2*u(M,N)...
        +d1*u(M+1,N+1)+d2*u(M,N+2);
elseif step==2
    a = -d2*ones((M-1)*N,1);
    b = (1+2*d2)*ones((M-1)*N,1);
    c = -d2*ones((M-1)*N,1);
    d = zeros((M-1)*N,1);
    % Case i=2 and j=2, includes BCs
    a(1)=0;
    b(1)=1+3*d2;
    d(1)=d1*u(3,2)+(1-2*d1)*u(2,2)+d1*u(1,2);
    %Case i=2 and j between [3,N], includes BCs
    for j=3:N
        d(j-1)=d1*u(3,j)+(1-2*d1)*u(2,j)+d1*u(1,j);
    end
    % Case i=2 and j=N+1, includes BCs
    b(N)=1+3*d2;
    c(N)=0;
    d(N)=d1*u(3,N+1)+(1-2*d1)*u(2,N+1)+d1*u(1,N+1);
    % Case i between [3,M-1] and j=2, includes BCs
    for i=3:M-1
        a((i-2)*N+1)=0;
        b((i-2)*N+1)=1+3*d2;
        d((i-2)*N+1)=d1*u(i+1,2)+(1-2*d1)*u(i,2)+d1*u(i-1,2);
    end
    % Interior of the interior, i in [3,M-1] and j in [3,N], no BCs
    for i=3:M-1
        for j=3:N
            d((i-2)*N+j-1)=d1*u(i+1,j)+(1-2*d1)*u(i,j)+d1*u(i-1,j);
        end
    end
    % Case i between [3,M-1] and j=N+1, includes BCs
    for i=3:M-1
        b((i-1)*N)=1+3*d2;
        c((i-1)*N)=0;
        d((i-1)*N)=d1*u(i+1,N+1)+(1-2*d1)*u(i,N+1)+d1*u(i-1,N+1);
    end
    % Case i=M and j=2, includes BCs
    a((M-2)*N+1)=0;
    b((M-2)*N+1)=1+3*d2;
    d((M-2)*N+1)=(1-2*d1)*u(M,2)+d1*u(M-1,2)...
        +d1*u(M+1,2);
    % Case i=M and j [3,N], includes BCs
    for j=3:N
        d((M-2)*N+j-1)=(1-2*d1)*u(M,j)+d1*u(M-1,j)...
            +d1*u(M+1,j);
    end
    % Case i=M and j=N+1, includes BCs
    b((M-1)*N)=1+3*d2;
    c((M-1)*N)=0;
    d((M-1)*N)=(1-2*d1)*u(M,N+1)+d1*u(M-1,N+1)...
        +d1*u(M+1,N+1);
end
end