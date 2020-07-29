function x=Utrisol(U,b)
% Column oriented version of upper triangular solver
n=length(b);
x=zeros(n,1);
for j=n:-1:2
    x(j)=b(j)/U(j,j);  %block multiplication
    b(1:j-1)=b(1:j-1)-x(j)*U(1:j-1,j);  %update b
end
x(1)=b(1)/U(1,1);
end