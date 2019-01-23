function [L,U,piv]=GEpiv(A)
[n,m]=size(A);
%First check that the matrix is square
if n~=m
    error('ERROR: The matrix given is not square.')
end
piv=(1:n)';
for k=1:n-1
    [maxV,r]=max(abs(A(k:n,k)));
    %If r=1 there is no row interchange.
    q=r+k-1; %position of the max in A(1:n,k)
    %Now interchange rows in the pivot vector and in the matrix A
    piv([k,q])=piv([q,k]); 
    A([k,q],:)=A([q,k],:);
    if A(k,k)~=0
        A(k+1:n,k)=A(k+1:n,k)/A(k,k);
        A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n); %Update A
    end
end
L=eye(n)+tril(A,-1);
U=triu(A);
end
        
    