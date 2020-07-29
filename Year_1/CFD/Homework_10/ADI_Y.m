function Y=ADI_Y(Y,M,N,dt,hx,hy,Re,Sc,Q)
if nargin<9
    disp('Source term zero')
    Q=zeros(size(Y));
end
%%%%%%%%%% STEP 1 %%%%%%%%%%
[a,b,c,d]=Yvectors(Y,M,N,dt,hx,hy,1/(Re*Sc),1,Q); %Obtain tridiagonal vectors
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% Correspondance with Y, obtain Y(n+1/2)
for j=4:N+3
    for i=4:M+3
        Y(i,j)=d((j-4)*M+i-3);
    end
end
% Boundary Conditions
Y=applyBCs(Y,M,N,hx,hy,'Y');

%%%%%%%%%% STEP 2 %%%%%%%%%%
[a,b,c,d]=Yvectors(Y,M,N,dt,hx,hy,1/(Re*Sc),2,Q); %Obtain tridiagonal vectors
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% Correspondance with Y, obtain Y(n+1)
for i=4:M+3
    for j=4:N+3
        Y(i,j)=d((i-4)*N+j-3);
    end
end
% Boundary Conditions
Y=applyBCs(Y,M,N,hx,hy,'Y');
end