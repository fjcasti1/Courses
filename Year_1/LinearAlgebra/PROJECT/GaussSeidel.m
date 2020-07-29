function [u,n]=GaussSeidel(N,tol)
% This function uses the Gauss-Seidel iteration method to come with the solution of the problem for given values of N and tol
h=1/N;          % Define the step
x=h*(0:N)';     % Create vector x
y=h*(0:N)';     % Create vector yuprev=rand(N+1,N+1);
% Create and initialize the solution to zero since is the mean value of the boundary conditions
u=zeros(N+1,N+1);
% Boundary conditions
for j=1:N+1
    u(1,j)=cos(pi*y(j));
    u(end,j)=exp(pi)*cos(pi*y(j));
    u(j,1)=exp(pi*x(j));
    u(j,end)=-exp(pi*x(j));
end
% Create the random matrix uprev to enter the while loop
uprev=rand(N+1,N+1);
% Initialize the counter of iterations n
n=0;
% While loop to keep iterating until the tolerance is met
while norm(uprev-u)>tol
    % Increase the value of n
    n=n+1;
    % Save the previous solution
    uprev=u;
    for i=2:N
        for j=2:N
            % Update the value of the matrix u with the values of the previous iteration and the updated values available
            u(i,j)=(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1))/4;
        end
    end
end