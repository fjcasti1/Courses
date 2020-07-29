function [fp_BFD,fp_CFD]=FiniteDifferences(N,f,x,h,gp)
    % This code gives the first order BFD and second order
    % CFD approximations for the function f. 
    % Arguments: -- N : size of the problem.
    %            -- f : function to approximate the derivative.
    %            -- x : vector with the nodes.   
    %            -- h : step.
    %            -- gp: if the mesh is not equispaced, and it is expressed by an analytical function.
    fp_BFD=zeros(N+1,1);
    fp_CFD=zeros(N+1,1);
    if nargin == 4
    % Backward Finite Differences
        fp_BFD(1)=(f(x(2))-f(x(1)))/h;
        for i=2:N+1
           fp_BFD(i)=(f(x(i))-f(x(i-1)))/h;
        end
    % Central Finite Differences
        fp_CFD(1)=(f(x(2))-f(x(1)))/h;
        for i=2:N
            fp_CFD(i)=(f(x(i+1))-f(x(i-1)))/(2*h);
        end
        fp_CFD(N+1)=(f(x(N+1))-f(x(N)))/h;
        
    elseif nargin == 5
        % Backward Finite Differences
        fp_BFD(1)=gp(x(1))*(f(x(2))-f(x(1)))/h;
        for i=2:N+1
           fp_BFD(i)=gp(x(i))*(f(x(i))-f(x(i-1)))/h;
        end
        % Central Finite Differences
        fp_CFD(1)=gp(x(1))*(f(x(2))-f(x(1)))/h;
        for i=2:N
            fp_CFD(i)=gp(x(i))*(f(x(i+1))-f(x(i-1)))/(2*h);
        end
        fp_CFD(N+1)=gp(x(N+1))*(f(x(N+1))-f(x(N)))/h;
    else 
        error('Wrong number of arguments')
    end
end