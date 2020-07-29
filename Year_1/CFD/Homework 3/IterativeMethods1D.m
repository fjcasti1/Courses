function [PHI,res]=IterativeMethods1D(f,phi_0,x,h,N,opt,iter,iter_res)
% This function calculates the solution of the 1D PDE given the initial
% guess phi_0, the right hand side f and information about the mesh.
% We use the parameters iter and iter_res to know what values to store and
% give back to the main code.
kmax=max([iter iter_res]);
j=1;
phi=phi_0(x);
PHI=zeros(length(phi),length(iter));
res=zeros(length(phi),iter_res+1);
if strcmp(opt,'PJ')
    % This includes the boundary conditions
    phi_new=zeros(size(phi));
    for k=1:kmax
        for i=2:N %Only the middle points will change value
            phi_new(i)=0.5*(phi(i+1)+phi(i-1)-h^2*f(x(i)));
        end
        phi=phi_new;
        if k == iter(j)
            PHI(:,j)=phi;
            j=j+1;
        end
        if k<=iter_res
            for i=2:N
                res(i,k+1)=f(x(i))-((phi(i+1)-2*phi(i)+phi(i-1))/h^2);
            end
        end
    end
elseif strcmp(opt,'GS')
    phi(1)=0;
    phi(end)=0;
    for k=1:kmax
        for i=2:N %Only the middle points will change value
            phi(i)=0.5*(phi(i+1)+phi(i-1)-h^2*f(x(i)));
        end
        if k == iter(j)
            PHI(:,j)=phi;
            j=j+1;
        end
        if k<=iter_res
            for i=2:N
                res(i,k+1)=f(x(i))-((phi(i+1)-2*phi(i)+phi(i-1))/h^2);
            end
        end
    end
elseif strcmp(opt,'SOR')
    phi(1)=0;
    phi(end)=0;
    w=2/(1+sqrt(1-(cos(pi/N))^2));
    for k=1:1:kmax
        for i=2:N %Only the middle points will change value
            phi(i)=phi(i)+w*(0.5*(phi(i+1)+phi(i-1)-h^2*f(x(i)))-phi(i));
        end
        if k == iter(j)
            PHI(:,j)=phi;
            j=j+1;
        end
        if k<=iter_res
            for i=2:N
                res(i,k+1)=f(x(i))-(phi(i+1)-2*phi(i)+phi(i-1))/h^2;
            end
        end
    end
end
for i=2:N
    res(i,1)=f(x(i))-((phi_0(x(i+1))-2*phi_0(x(i))+phi_0(x(i-1)))/h^2);
end
end