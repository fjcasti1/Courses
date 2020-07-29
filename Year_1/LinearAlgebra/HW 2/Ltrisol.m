function x=Ltrisol(L,b)
%Column orientated version of lower diagonal solver
n=length(b);
x=zeros(n,1); %reserved storage of x
for j=1:n-1
    x(j)=b(j)/L(j,j); %block multiplication
    b(j+1:n)=b(j+1:n)-x(j)*L(j+1:n,j);  %update b
end
x(n)=b(n)/L(n,n);
end