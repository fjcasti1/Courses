function [x, niter, conv] = newton_vector(f,df,x0,tol,maxiter)
% Newton's method for rootfinding. Multivariable case.
%   NEWTON(F,DF,X0,TOL,MAXITER) returns the root of the function handle F
%   DF(x) returns the Jacobian of F at a point x
%   X0 is an initial guess
%   TOL is a tolerance for breaking the loop
%   MAXITER is the maximum number of iterations allowed
%
%   Output: [X, NITER, CONV]
%   X is the solution 
%   NITER is the number of iterations (optional)
%   CONV is true if convergence was successful and false otherwise
%   (optional)
%
%   Rodrigo Platte, 2008

%% Initialization of variables
x(1) = x0;             % initial guess   
dx = inf;           % difference between iterates, set INF to force loop
niter = 0;          % number of iterations
conv = true;        % assume convergence
if nargin<5         % if maxiter isn't provided as input
    maxiter = 2e3;  % maximum number of iterations allowed
end

%% Main loop
while norm(dx,inf)>tol %&& norm(f(x),inf)>tol

    % loop breaks if number of iterations exceeds maxiter
    niter = niter+1;
    if niter == maxiter
        conv = false;           % failed to converge
%        warning('Newton1:niter','Reached maximum number of iterations')
        break
    end
    
    %xold = x(niter);                   % update previous iterate
    
    x(niter+1) = x(niter) - df(x(niter))\f(x(niter));         % Newton step
    dx = x(niter+1) - x(niter);              % difference between iterates  
end

% Display number of steps used to obtain the answer
%if niter~=maxiter
%    disp(['converged in ', int2str(niter-1),' steps'])
%end