function u=ADI_u(u,M,N,dt,time,Q)
%%%%%%%%%% STEP 1 %%%%%%%%%%
[a,b,c,d]=uvectors(u,M,N,dt,1,Q); %Obtain tridiagonal vectors
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% Correspondance with u, obtain u(n+1/2)
for j=2:N+1
    for i=2:M
        u(i,j)=d((j-2)*(M-1)+i-1);
    end
end
% Update top and bottom boundaries (ghost cells), they depend on the
% interior, I don't call the BCs function because the for loops in it are
% unnecessary.
% u(:,N+2)=-u(:,N+1);
% u(:,1)=-u(:,2);
u=applyBCs(u,M,N,time,'u');

%%%%%%%%%% STEP 2 %%%%%%%%%%
[a,b,c,d]=uvectors(u,M,N,dt,2,Q); %Obtain tridiagonal vectors
d=GaussTriSol(a,b,c,d);
for i=2:M
    for j=2:N+1
        u(i,j)=d((i-2)*N+j-1);
    end
end
% Update top and bottom boundaries (ghost cells), they depend on the
% interior, I don't call the BCs function because the for loops in it are
% unnecessary.
% u(:,N+2)=-u(:,N+1); % Top
% u(:,1)=-u(:,2); % Bottom
u=applyBCs(u,M,N,time,'u');
end