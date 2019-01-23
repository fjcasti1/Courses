%% Test for Householder Decomposition
clear all
close all
clc
%Define a random matrix and check that is non-singular with the condition
%number of A
m=5;
n=5;
A=zeros(m,n);
tol=100;
i=1;
%The following loop will redefine the matrix A until one of them has a
%condition number below tol.
while cond(A)>tol
A=rand(m,n); %Redefine A
i=i+1;      %Count how many tries
end
% A=[0 0 3 2; 1 2 1 7; 0 0 0 1; 0 2 1 3; 1 1 1 1];
% A=[0 0 3; 1 2 1; 0 0 0; 0 2 1];
% A=[1 2 3; 4 5 6; 7 8 9];
% A=[1 2 3 4 5; 6 7 8 9 10; 95 94 93 92 91; 55 56 57 58 59; 30 40 50 60 70];
m=size(A,1)
n=size(A,2)
A_0=A;
if m==n
    nfinal=n-1;
elseif m>n
    nfinal=n;
end
for j=1:nfinal
[v,beta]=house(A(j:m,j));
A(j:m,j:n)=A(j:m,j:n)-(beta*v)*(v'*A(j:m,j:n));
A(j+1:m,j)=v(2:m-j+1);
end
Q=eye(m);
for j=nfinal:-1:1
    v(j:m)=[1;A(j+1:m,j)];
    beta=2/(1+norm(A(j+1:m,j))^2);
    Q(j:m,j:m)=Q(j:m,j:m)-(beta*v(j:m))*(v(j:m)'*Q(j:m,j:m));
end
% end
R=triu(A);
[q,r]=qr(A_0);
q
Q
r
R
q*r
Q*R
A_0