function [A,P,D]=DominantEigenvalueMatrix(N,f)
% This function gives a NxN matrix A, an orthogonal matrix P and a diagonal
% matrix D such that A=A = P*D*P'. Needs the dimension N and a factor f. 
% The matrix A will have an dominant eigenvalue if f>>1 and will have the two
% larger eigenvalues very similar if f is close to unity.
P = orth(rand(N));
lambdaV = randi([1,100],N,1);
k=randi([1,N],1);
j=find(lambdaV==max(lambdaV));
while k==j
    k=randi([1,N],1);
end
lambdaV(k)=f*max(lambdaV);
D = diag(lambdaV);
A = P*D*P';
end