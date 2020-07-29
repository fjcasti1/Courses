function [a,b,c,d] = vvectors(v,M,N,dt,hx,hy,alpha,step)
xi=linspace(-hx/2,4+hx/2,M+2);
yj=linspace(0,2,N+1);
dx=alpha*dt/hx^2;
dy=alpha*dt/hy^2;
d1=dx/2; %debug d1=1/8;
d2=dy/2;
if step==1
    a = -d1*ones(M*(N-1),1);
    b = (1+2*d1)*ones(M*(N-1),1);
    c = -d1*ones(M*(N-1),1);
    d = zeros(M*(N-1),1);
    % Case j=2 and i=2, includes BCs but they are already
    % imposed for the current timestep
    a(1)=0;
    b(1)=1+3*d1;
    d(1)=d2*v(2,3)+(1-2*d2)*v(2,2)+d2*v(2,1);
    % Case j=2 and i between [3,M], includes BCs but they are already
    % imposed for the current timestep
    for i=3:M
        d(i-1)=d2*v(i,3)+(1-2*d2)*v(i,2)+d2*v(i,1); 
    end
    % Case j=2 and i=M+1, includes BCs
    b(M)=1+3*d1;
    c(M)=0;
    d(M)=d2*v(M+1,3)+(1-2*d2)*v(M+1,2)+d2*v(M+1,1);
    % Case j between [3,N-1] and i=2, includes BCs
    for j=3:N-1
        a((j-2)*M+1)=0;
        b((j-2)*M+1)=1+3*d1;
        d((j-2)*M+1)=d2*v(2,j+1)+(1-2*d2)*v(2,j)+d2*v(2,j-1);
    end
    % Interior of the interior, j in [3,N-1] and i in [3,M], no BCs
    for j=3:N-1
        for i=3:M
            d((j-2)*M+i-1)=d2*v(i,j+1)+(1-2*d2)*v(i,j)+d2*v(i,j-1);
        end
    end
    % Case j between [3,N-1] and i=M+1, includes BCs
    for j=3:N-1
        b((j-1)*M)=1+3*d1;
        c((j-1)*M)=0;
        d((j-1)*M)=d2*v(M+1,j+1)+(1-2*d2)*v(M+1,j)+d2*v(M+1,j-1);
    end
    % Case j=N and i=2, includes BCs but they are already
    % imposed for the current timestep
    a((N-2)*M+1)=0;
    b((N-2)*M+1)=1+3*d1;
    d((N-2)*M+1)=d2*v(2,N+1)+(1-2*d2)*v(2,N)+d2*v(2,N-1);
    % Case j=N and i [3,M], includes BCs but they are already
    % imposed for the current timestep
    for i=3:M
        d((N-2)*M+i-1)=d2*v(i,N+1)+(1-2*d2)*v(i,N)+d2*v(i,N-1);
    end
    % Case j=N and i=M+1, includes BCs but they are already
    % imposed for the current timestep
    b((N-1)*M)=1+3*d1;
    c((N-1)*M)=0;
    d((N-1)*M)=d2*v(M+1,N+1)+(1-2*d2)*v(M+1,N)+d2*v(M+1,N-1);
elseif step==2
    a = -d2*ones(M*(N-1),1);
    b = (1+2*d2)*ones(M*(N-1),1);
    c = -d2*ones(M*(N-1),1);
    d = zeros(M*(N-1),1);
    % Case j=2 and i=2, includes BCs
    a(1)=0;
    d(1)=d1*v(3,2)+(1-3*d1)*v(2,2);
    %Case i=2 and j between [3,N-1], includes BCs
    for j=3:N-1
        d(j-1)=d1*v(3,j)+(1-3*d1)*v(2,j);
    end
    % Case i=2 and j=N, includes BCs
    c(N-1)=0;
    d(N-1)=d1*v(3,N)+(1-3*d1)*v(2,N)+d2*v(2,N+1);
    % Case i between [3,M] and j=2, includes BCs
    for i=3:M
        a((i-2)*(N-1)+1)=0;
        if (xi(i)>=1.5 && xi(i)<=2.5)
            b((i-2)*(N-1)+1)=1+2*d2/3;
            c((i-2)*(N-1)+1)=-2*d2/3;
        end
        d((i-2)*(N-1)+1)=d1*v(i+1,2)+(1-2*d1)*v(i,2)+d1*v(i-1,2);
    end
    % Interior of the interior, i in [3,M] and j in [3,N-1], no BCs
    for i=3:M
        for j=3:N-1
            d((i-2)*(N-1)+j-1)=d1*v(i+1,j)+(1-2*d1)*v(i,j)+d1*v(i-1,j);
        end
    end
    % Case i between [3,M] and j=N, includes BCs
    for i=3:M
        c((i-1)*(N-1))=0;
        d((i-1)*(N-1))=d1*v(i+1,N)+(1-2*d1)*v(i,N)+d1*v(i-1,N)...
            +d2*v(i,N+1);
    end
    % Case i=M+1 and j=2, includes BCs
    a((M-1)*(N-1)+1)=0;
    d((M-1)*(N-1)+1)=(1-3*d1)*v(M+1,2)+d1*v(M,2);
    % Case i=M+1 and j [3,N-1], includes BCs
    for j=3:N-1
        d((M-1)*(N-1)+j-1)=(1-3*d1)*v(M+1,j)+d1*v(M,j);
    end
    % Case i=M+1 and j=N, includes BCs
    c((N-1)*M)=0;
    d((N-1)*M)=(1-3*d1)*v(M+1,N)+d1*v(M,N)...
        +d2*v(M+1,N+1);
end
end