function [lambda,k,q]=powermethod(A,tol)
% This function uses the powermethod to, given the matrix A and a tolerance,
% obtain the eigenvalue with larger absolute value and its eigenvector. 
% It will also provide the number of iterations needed to meet the tolerance.
    N=size(A,1);
    lambdaprev=1;   % Initialize lambdaprev
    lambda=0;       % Initialize lambda
    k=0;            % Start the counter of iterations
    q=rand(N,1);    % The first guess of q is a random vector as the problem specifies
    while norm(lambdaprev-lambda)>tol   % This is the power method algorithm to obtain the dominant eigenvalue
        k=k+1;
        lambdaprev=lambda;
        z=A*q;
        q=z/norm(z);
        lambda=q'*A*q;
    end
end