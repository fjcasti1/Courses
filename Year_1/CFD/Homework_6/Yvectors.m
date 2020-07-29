function [a,b,c,d] = Yvectors(Y,M,N,dt,hx,hy,alpha,step)
xi=linspace(-hx/2,4+hx/2,M+2);
yj=linspace(-hy/2,2+hy/2,N+2);
dx=alpha*dt/hx^2;
dy=alpha*dt/hy^2;
d1=dx/2;
d2=dy/2;

if step==1
    a = -d1*ones(M*N,1);
    b = (1+2*d1)*ones(M*N,1);
    c = -d1*ones(M*N,1);
    d = zeros(M*N,1);
    % Case j=2 and i=2, includes BCs but they are already
    % imposed for the current timestep
    a(1)=0;
    b(1)=1+d1;
    d(1)=d2*Y(2,3)+(1-2*d2)*Y(2,2)+d2*Y(2,1);
    % Case j=2 and i between [3,M], includes BCs but they are already
    % imposed for the current timestep
    for i=3:M
        d(i-1)=d2*Y(i,3)+(1-2*d2)*Y(i,2)+d2*Y(i,1); 
    end
    % Case j=2 and i=M+1, includes BCs
    b(M)=1+d1;
    c(M)=0;
    d(M)=d2*Y(M+1,3)+(1-2*d2)*Y(M+1,2)+d2*Y(M+1,1);
    % Case j between [3,N] and i=2, includes BCs 
    for j=3:N
        a((j-2)*M+1)=0;
        if (yj(j)>=0.5 && yj(j)<=1)
            b((j-2)*M+1)=1+3*d1;            
            d((j-2)*M+1)=d2*Y(2,j+1)+(1-2*d2)*Y(2,j)+d2*Y(2,j-1)...
                +2*d1;
        else
            b((j-2)*M+1)=1+d1;
            d((j-2)*M+1)=d2*Y(2,j+1)+(1-2*d2)*Y(2,j)+d2*Y(2,j-1);
        end
    end
    % Interior of the interior, j in [3,N] and i in [3,M], no BCs
    for j=3:N
        for i=3:M
            d((j-2)*M+i-1)=d2*Y(i,j+1)+(1-2*d2)*Y(i,j)+d2*Y(i,j-1);
        end
    end
    % Case j between [3,N] and i=M+1, includes BCs
    for j=3:N
        c((j-1)*M)=0;
        if (yj(j)>=1 && yj(j)<=1.5)
            b((j-1)*M)=1+3*d1;            
            d((j-1)*M)=d2*Y(M+1,j+1)+(1-2*d2)*Y(M+1,j)+d2*Y(M+1,j-1)...
                +0.5*d1;
        else
            b((j-1)*M)=1+d1;
            d((j-1)*M)=d2*Y(M+1,j+1)+(1-2*d2)*Y(M+1,j)+d2*Y(M+1,j-1);
        end
    end
    % Case j=N+1 and i=2, includes BCs but they are already
    % imposed for the current timestep
    a((N-1)*M+1)=0;
    b((N-1)*M+1)=1+d1;
    d((N-1)*M+1)=d2*Y(2,N+2)+(1-2*d2)*Y(2,N+1)+d2*Y(2,N);
    % Case j=N+1 and i [3,M], includes BCs but they are already
    % imposed for the current timestep
    for i=3:M
        d((N-1)*M+i-1)=d2*Y(i,N+2)+(1-2*d2)*Y(i,N+1)+d2*Y(i,N);
    end
    % Case j=N+1 and i=M+1, includes BCs but they are already
    % imposed for the current timestep
    b(N*M)=1+d1;
    c(N*M)=0;
    d(N*M)=d2*Y(M+1,N+2)+(1-2*d2)*Y(M+1,N+1)+d2*Y(M+1,N);
elseif step==2
    a = -d2*ones(M*N,1);
    b = (1+2*d2)*ones(M*N,1);
    c = -d2*ones(M*N,1);
    d = zeros(M*N,1);
    % Case j=2 and i=2, includes BCs
    a(1)=0;
    b(1)=1+d2;
    d(1)=d1*Y(3,2)+(1-2*d1)*Y(2,2)+d1*Y(1,2);
    %Case i=2 and j between [3,N], includes BCs
    for j=3:N
        d(j-1)=d1*Y(3,j)+(1-2*d1)*Y(2,j)+d1*Y(1,j);
    end
    % Case i=2 and j=N+1, includes BCs
    b(N)=1+d2;
    c(N)=0;
    d(N)=d1*Y(3,N+1)+(1-2*d1)*Y(2,N+1)+d1*Y(1,N+1);
    % Case i between [3,M] and j=2, includes BCs
    for i=3:M
        a((i-2)*N+1)=0;
        b((i-2)*N+1)=1+d2;
        d((i-2)*N+1)=d1*Y(i+1,2)+(1-2*d1)*Y(i,2)+d1*Y(i-1,2);
    end
    % Interior of the interior, i in [3,M] and j in [3,N], no BCs
    for i=3:M
        for j=3:N
            d((i-2)*N+j-1)=d1*Y(i+1,j)+(1-2*d1)*Y(i,j)+d1*Y(i-1,j);
        end
    end
    % Case i between [3,M] and j=N+1, includes BCs
    for i=3:M
        c((i-1)*N)=0;
        if (xi(i)>=0.5 && xi(i)<=1)
            b((i-1)*N)=1+3*d2;
        else
            b((i-1)*N)=1+d2;
        end
        d((i-1)*N)=d1*Y(i+1,N+1)+(1-2*d1)*Y(i,N+1)+d1*Y(i-1,N+1);
    end
    % Case i=M+1 and j=2, includes BCs
    a((M-1)*N+1)=0;
    b((M-1)*N+1)=1+d2;
    d((M-1)*N+1)=d1*Y(M+2,2)+(1-2*d1)*Y(M+1,2)+d1*Y(M,2);
    % Case i=M+1 and j [3,N], includes BCs
    for j=3:N
        d((M-1)*N+j-1)=d1*Y(M+2,j)+(1-2*d1)*Y(M+1,j)+d1*Y(M,j);
    end
    % Case i=M+1 and j=N+1, includes BCs
    b(M*N)=1+d2;
    c(M*N)=0;
    d(M*N)=d1*Y(M+2,N+1)+(1-2*d1)*Y(M+1,N+1)+d1*Y(M,N+1);
end
end