function Y=ADI_Y(Y,M,N,dt,Q)
%%%%%%%%%% STEP 1 %%%%%%%%%%
[a,b,c,d]=Yvectors(Y,M,N,dt,1,Q); %Obtain tridiagonal vectors
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% Correspondance with Y, obtain Y(n+1/2)
for j=4:N+3
    for i=4:M+3
        Y(i,j)=d((j-4)*M+i-3);
    end
end
% Boundary Conditions
Y=applyBCs(Y,M,N,'Y');

%%%%%%%%%% STEP 2 %%%%%%%%%%
[a,b,c,d]=Yvectors(Y,M,N,dt,2,Q); %Obtain tridiagonal vectors
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% Correspondance with Y, obtain Y(n+1)
for i=4:M+3
    for j=4:N+3
        Y(i,j)=d((i-4)*N+j-3);
    end
end
% Boundary Conditions
Y=applyBCs(Y,M,N,'Y');
end