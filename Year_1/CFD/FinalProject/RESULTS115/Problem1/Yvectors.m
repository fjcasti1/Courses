function [a,b,c,d] = Yvectors(Y,M,N,dt,step,Q)
global Re;
global Sc;
global hx;
global hy;
alpha=1/(Re*Sc);

Yinlet1=1;
Yinlet2=1;

dx=alpha*dt/hx^2;
dy=alpha*dt/hy^2;
d1=dx/2;
d2=dy/2;
x=linspace(-5*hx/2,4+5*hx/2,M+6);
y=linspace(-5*hy/2,2+5*hy/2,N+6);
if step==1
    a = -d1*ones(M*N,1);
    b = (1+2*d1)*ones(M*N,1);
    c = -d1*ones(M*N,1);
    d = zeros(M*N,1);
    for j=4:N+3
        for i=4:M+3
            if i==4 % Left boundary BCs
                a((j-4)*M+1)=0; 
                if (y(j)>=0.5 && y(j)<=1)
                    b((j-4)*M+1)=1+3*d1; 
                    d((j-4)*M+1)=d2*Y(4,j-1)+(1-2*d2)*Y(4,j)+d2*Y(4,j+1)...
                                    +d1*2*Yinlet1+(dt/2)*Q(i,j); 
                else
                    b((j-4)*M+1)=1+d1; 
                    d((j-4)*M+1)=d2*Y(4,j-1)+(1-2*d2)*Y(4,j)+d2*Y(4,j+1)...
                        +(dt/2)*Q(i,j);
                end
            elseif i==M+3% Right boundary BCs
                c((j-4)*M+M)=0;
                if (y(j)>=1 && y(j)<=1.5)
                    b((j-4)*M+M)=1+3*d1;
                    d((j-4)*M+M)=d2*Y(i,j-1)+(1-2*d2)*Y(i,j)+d2*Y(i,j+1)...
                                    +d1*2*Yinlet2+(dt/2)*Q(i,j); 
                else
                    b((j-4)*M+M)=1+d1; 
                    d((j-4)*M+M)=d2*Y(i,j-1)+(1-2*d2)*Y(i,j)+d2*Y(i,j+1)...
                                    +(dt/2)*Q(i,j);
                end
            else % The top and bottom boundaries are imposed by updating the ghost
                    % cells accordingly with them.
                d((j-4)*M+i-3)=d2*Y(i,j-1)+(1-2*d2)*Y(i,j)+d2*Y(i,j+1)+(dt/2)*Q(i,j);
            end
        end
    end
elseif step==2
    a = -d2*ones(M*N,1);
    b = (1+2*d2)*ones(M*N,1);
    c = -d2*ones(M*N,1);
    d = zeros(M*N,1);
    for i=4:M+3
        for j=4:N+3
            if j==4 % Bottom boundary BCs
                a((i-4)*N+1)=0;
                b((i-4)*N+1)=1+d2;
                d((i-4)*N+1)=d1*Y(i-1,4)+(1-2*d1)*Y(i,4)+d1*Y(i+1,4)...
                                +(dt/2)*Q(i,j);
            elseif j==N+3 % Top boundary BCs
                c((i-4)*N+N)=0;
                if (x(i)>=0.5 && x(i)<=1)
                    b((i-4)*N+N)=1+3*d2;
                else
                    b((i-4)*N+N)=1+d2;
                end
                d((i-4)*N+N)=d1*Y(i-1,N+3)+(1-2*d1)*Y(i,N+3)+d1*Y(i+1,N+3)...
                                +(dt/2)*Q(i,j);
            else % The left and right boundaries are imposed by initialization,
                    % they do not change
                d((i-4)*N+j-3)=d1*Y(i-1,j)+(1-2*d1)*Y(i,j)+d1*Y(i+1,j)...
                                +(dt/2)*Q(i,j);
            end
        end
    end
end
end