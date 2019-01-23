function v=ADI_v(v,M,N,dt,hx,hy,Re,Q)
if nargin<8
    disp('Source term zero')
    Q=zeros(size(u));
end
%%%%%%%%%% STEP 1 %%%%%%%%%%
[a,b,c,d]=vvectors(v,M,N,dt,hx,hy,1/Re,1,Q); %Obtain tridiagonal vectors
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% Correspondance with v, obtain v(n+1/2)
for j=2:N
    for i=2:M+1
        v(i,j)=d((j-2)*M+i-1);
    end
end
% Update left and right boundaries (ghost cells), as well as the outlet
% Neumann BC, they depend on the interior.
v=applyBCs(v,M,N,hx,hy,'v');

%%%%%%%%%% STEP 2 %%%%%%%%%%
[a,b,c,d]=vvectors(v,M,N,dt,hx,hy,1/Re,2,Q); %Obtain tridiagonal vectors
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% Correspondance with v, obtain v(n+1/2)
for i=2:M+1
    for j=2:N
        v(i,j)=d((i-2)*(N-1)+j-1);
    end
end
% Update left and right boundaries (ghost cells), as well as the outlet
% Neumann BC, they depend on the interior.
v=applyBCs(v,M,N,hx,hy,'v');
end