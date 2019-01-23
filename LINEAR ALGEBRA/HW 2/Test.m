%% Test for Gaussian elimination and solvers

%Define a random matrix and check that is non-singular with the condition
%number of A
A=zeros(20,20);
tol=100;
i=1;
%The following loop will redefine the matrix A until one of them has a
%condition number below tol.
while cond(A)>tol
A=rand(20); %Redefine A
i=i+1;      %Count how many tries
end
A
i
b=rand(20,1)
%We proceed with the Gaussian Elimination with Partial Pivoting.
[L,U,piv]=GEpiv(A)
%Then we obtain PA=A(piv,:), L and U, where PA=LU
A(piv,:)
L*U
%Observe that they are equal
%We are going to solve the following system: PAx=LUx=Pb, where Pb=b(piv).
%We start with the following: Lz=Pb, z=Ux.
z=Ltrisol(L,b(piv))
%Obtained z, we now calculate the solution x
x=Utrisol(U,z)
%Calculate the residue vector
r=b(piv)-A(piv,:)*x
%Calculate two norms of the residue vector
N1=norm(r,1)
N2=norm(r,2)

% About the row interchanges, the information is given by the vector piv.