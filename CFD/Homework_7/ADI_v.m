function v=ADI_v(v,M,N,dt,hx,hy,Re)
xiv=linspace(-hx/2,4+hx/2,M+2);
%%%%%%%%%% STEP 1 %%%%%%%%%%
[a,b,c,d]=vvectors(v,M,N,dt,hx,hy,1/Re,1); %Obtain tridiagonal vectors
% keyboard
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% keyboard
% Correspondance with v, obtain v(n+1/2)
for j=2:N
    for i=2:M+1
        v(i,j)=d((j-2)*M+i-1);
    end
end
% Update left and right boundaries (ghost cells), as well as the outlet
% Neumann BC, they depend on the interior.
v(1,:)=-v(2,:); % Left
v(M+2,:)=-v(M+1,:); % Right
for i=1:M+2
    if (xiv(i)>=1.5 && xiv(i)<=2.5) % Outlet
        v(i,1)=(4/3)*v(i,2)-(1/3)*v(i,3);
    end
end

%%%%%%%%%% STEP 2 %%%%%%%%%%
[a,b,c,d]=vvectors(v,M,N,dt,hx,hy,1/Re,2); %Obtain tridiagonal vectors
% keyboard
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% keyboard
% Correspondance with v, obtain v(n+1/2)
for i=2:M+1
    for j=2:N
        v(i,j)=d((i-2)*(N-1)+j-1);
    end
end
% Update left and right boundaries (ghost cells), as well as the outlet
% Neumann BC, they depend on the interior.
v(1,:)=-v(2,:); % Left
v(M+2,:)=-v(M+1,:); % Right
for i=1:M+2
    if (xiv(i)>=1.5 && xiv(i)<=2.5) % Outlet
        v(i,1)=(4/3)*v(i,2)-(1/3)*v(i,3);
    end
end
end