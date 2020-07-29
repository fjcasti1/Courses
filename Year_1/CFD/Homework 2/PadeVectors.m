function [a,b,c,d]=PadeVectors(N,f,x,h)
    % This function gives four vectors corresponding to the three diagonals
    % of the PADE scheme and the column vector of coefficients.
    % Arguments: -- N : size of the problem.
    %            -- f : function to approximate the derivative.
    %            -- x : vector with the nodes.   
    %            -- h : step.
    
    a=ones(N+1,1);a(end)=2;
    b=4*ones(N+1,1);
    b(1)=1;b(end)=1;
    c=ones(N+1,1);c(1)=2;
    d=ones(N+1,1);
    d(1)=-2.5*f(x(1))+2*f(x(2))+0.5*f(x(3));
    for i=2:N
        d(i)=3*(f(x(i+1))-f(x(i-1)));
    end
    d(N+1)=2.5*f(x(N+1))-2*f(x(N))-0.5*f(x(N-1));    
    d=d/h;
end