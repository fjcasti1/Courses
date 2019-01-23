%% FRANCISCO CASTILLO APM 505 PROJECT

%% Initialization of the code
clear all    % Clear workspace
close all    % Close all figures
clc          % Clear command window
% Create four column vectors for storing the number of iterations for each method.
nGJ=zeros(5:1);
nGS=zeros(5:1);
nSOR=zeros(5:1);
normGJ=zeros(5:1);
normGS=zeros(5:1);
normSOR=zeros(5:1);
% Assign values of N and tol for each one of the four cases given in the assignment.
for i=1:5
    switch i
        case 1
            N=8;
            tol=1e-3;
        case 2
            N=16;
            tol=1e-4;
        case 3
            N=32;
            tol=1e-5;
        case 4
            N=64;
            tol=1e-6;
        case 5
            N=128;
            tol=1e-7;
    end
    
%% Obtain analytical solution
x=linspace(0,1,N+1);
[Y,X]=meshgrid(x,x);
U=exp(pi*X).*cos(pi*Y);

%% Run the different iteration methods    
[uGJ,nGJ(i)]=GaussJacobi(N,tol);
[uGS,nGS(i)]=GaussSeidel(N,tol);
[uSOR,nSOR(i)]=SOR(N,tol);

%% Store the difference with the analytical solution
normGJ(i)=norm(uGJ-U);
normGS(i)=norm(uGS-U);
normSOR(i)=norm(uSOR-U);
end
%% Present the results in a table
T=table;
T.N = [8;16;32;64;128];
T.tol = ['1e-3';'1e-4';'1e-5';'1e-6';'1e-7'];
T.nGJ= nGJ';
T.nGS= nGS';
T.nSOR= nSOR';
T.NormGJ=normGJ';
T.NormGS=normGS';
T.NormSOR=normSOR'

%% Function for the Gauss-Jacobi iteration method

type GaussJacobi.m

%% Function for the Gauss-Seidel iteration method

type GaussSeidel.m

%% Function for the SOR iteration method

type SOR.m