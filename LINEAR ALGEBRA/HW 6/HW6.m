%% FRANCISCO CASTILLO APM 505 HOMEWORK 6

%% Initialization of the code
clear all    % Clear workspace
clc          % Clear command window
format long

tol=1e-12; % Tolerance
%% Function to create the matrix A
    type DominantEigenvalueMatrix.m
%% Function for the Power Method Iteration
    type powermethod.m
%% Run different cases of study
for i=1:4
    switch i
        case 1
            N=3;    % Dimension of the matrix A
            f=30;   % The factor f large means that the matrix A is going to have one dominant eigenvalue
        case 2
            N=3;
            f=1.0001;   % The factor f close to unity means that the matrix A is not going to have any dominant eigenvalue
        case 3
            N=9;
            f=30;
        case 4
            N=9;
            f=1.0001;
    end
%% Create the matrix A
    [A,P,D]=DominantEigenvalueMatrix(N,f);
%% Power Method Iteration
    [lambda,k,q]=powermethod(A,tol);
%% Results and discussion
    fprintf('>>Case %d\n',i)
    A
    P
    D
    lambda
    k
    q
    v=P'*q
end