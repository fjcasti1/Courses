function [PHI,res]=IterativeMethods2D(f,phi_0,x,h,M,N,opt,iter,iter_res)
% This function calculates the solution of the 2D PDE given the initial
% guess phi_0, the right hand side f and information about the mesh.
% We use the parameters iter and iter_res to know what values to store and
% give back to the main code.
kmax=max([iter iter_res]);
l=1;
phi=zeros(N+1,N+1);
for i=1:N+1
    for j=1:N+1
        phi(i,j)=phi_0(x(i),x(j));
    end
end
PHI=zeros(length(phi),length(phi),length(iter));% Includes BCs
res=zeros(length(phi),length(phi),iter_res+1);
if strcmp(opt,'GS')
    for k=1:kmax
        for i=2:N %Only the middle points will change value
            for j=2:N
                phi(i,j)=0.25*(phi(i+1,j)+phi(i-1,j)...
                    +phi(i,j+1)+phi(i,j-1)-h^2*f(x(i),x(j)));
            end
        end
        if k == iter(l)           
            PHI(:,:,l)=phi;
            l=l+1;
        end
        if k<=iter_res
            for i=2:N
                for j=2:N
                    res(i,j,k+1)=f(x(i),x(j))-(phi(i+1,j)-2*phi(i,j)+phi(i-1,j))/h^2-...
                        (phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/h^2;
                end
            end
        end
    end
elseif strcmp(opt,'SOR')
    rho=0.5*(cos(pi/M)+cos(pi/N));
    w=2/(1+sqrt(1-rho^2));
    for k=1:kmax
        for i=2:N %Only the middle points will change value
            for j=2:N
                phi(i,j)=phi(i,j)+w*(0.25*(phi(i+1,j)+phi(i-1,j)...
                    +phi(i,j+1)+phi(i,j-1)-h^2*f(x(i),x(j)))-phi(i,j));
            end
        end
        if k == iter(l)           
            PHI(:,:,l)=phi;
            l=l+1;
        end
        if k<=iter_res
            for i=2:N
                for j=2:N
                    res(i,j,k+1)=f(x(i),x(j))-(phi(i+1,j)-2*phi(i,j)+phi(i-1,j))/h^2-...
                        (phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/h^2;
                end
            end
        end
    end
end
for i=2:N
    for j=2:N
        res(i,j,1)=f(x(i),x(j))-(phi_0(i+1,j)-2*phi_0(i,j)+phi_0(i-1,j))/h^2-...
                        (phi_0(i,j+1)-2*phi_0(i,j)+phi_0(i,j-1))/h^2;
    end
end
end