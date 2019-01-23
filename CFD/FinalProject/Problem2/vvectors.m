function [a,b,c,d] = vvectors(v,M,N,dt,step,Q)
global hx;
global hy;
global xv;
global Re;
alpha=1/Re;

global ao; global bo;

dx=alpha*dt/hx^2;
dy=alpha*dt/hy^2;
d1=dx/2;
d2=dy/2;
if step==1
    a = -d1*ones(M*(N-1),1);
    b = (1+2*d1)*ones(M*(N-1),1);
    c = -d1*ones(M*(N-1),1);
    d = zeros(M*(N-1),1);
    for j=2:N
        for i=2:M+1
            if i==2
                a((j-2)*M+1)=0;
                b((j-2)*M+1)=1+3*d1;
                d((j-2)*M+1)=d2*v(2,j-1)+(1-2*d2)*v(2,j)+d2*v(2,j+1)+(dt/2)*Q(i,j);
            elseif i==M+1
                b((j-1)*M)=1+3*d1;
                c((j-1)*M)=0;
                d((j-1)*M)=d2*v(M+1,j-1)+(1-2*d2)*v(M+1,j)+d2*v(M+1,j+1)+(dt/2)*Q(i,j);
            elseif (j==2 && xv(i)>=ao && xv(i)<=bo)
                d(i-1)=d2*(4*v(i,2)/3-v(i,3)/3)+(1-2*d2)*v(i,2)+d2*v(i,3)+(dt/2)*Q(i,j);
            else
                d((j-2)*M+i-1)=d2*v(i,j-1)+(1-2*d2)*v(i,j)+d2*v(i,j+1)+(dt/2)*Q(i,j);
            end
        end
    end
elseif step==2
    a = -d2*ones(M*(N-1),1);
    b = (1+2*d2)*ones(M*(N-1),1);
    c = -d2*ones(M*(N-1),1);
    d = zeros(M*(N-1),1);
    for i=2:M+1
        for j=2:N
            if j==2
                a((i-2)*(N-1)+1)=0;
                if (xv(i)>=ao && xv(i)<=bo)
                    b((i-2)*(N-1)+1)=1+(2-4/3)*d2;
                    c((i-2)*(N-1)+1)=-(2/3)*d2;
                end
                d((i-2)*(N-1)+1)=d1*v(i-1,2)+(1-2*d1)*v(i,2)+d1*v(i+1,2)+(dt/2)*Q(i,j);
            elseif j==N
                c((i-1)*(N-1))=0;
                d((i-1)*(N-1))=d1*v(i-1,N)+(1-2*d1)*v(i,N)+d1*v(i+1,N)...
                                +(dt/2)*Q(i,j)+d2*v(i,N+1);
            else
                d((i-2)*(N-1)+j-1)=d1*v(i-1,j)+(1-2*d1)*v(i,j)+d1*v(i+1,j)+(dt/2)*Q(i,j);
            end
        end
    end
end
end