function Y=ADI_Y(Y,M,N,dt,hx,hy,Re,Sc)
%% Solve for Y(x,y,t)
xiY=linspace(-hx/2,4+hx/2,M+2);
yjY=linspace(-hy/2,2+hy/2,N+2);
%%%%%%%%%% STEP 1 %%%%%%%%%%
[a,b,c,d]=Yvectors(Y,M,N,dt,hx,hy,1/(Re*Sc),1); %Obtain tridiagonal vectors
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% Correspondance with Y, obtain Y(n+1/2)
for j=2:N+1
    for i=2:M+1
        Y(i,j)=d((j-2)*M+i-1);
    end
end
% Update ghost cells.
Y(:,1)=Y(:,2); % Bottom
Y(:,N+2)=Y(:,N+1); % Top, we will overwrite the inlet 3
Y(find(xiY>=0.5,1):find(xiY<=1,1,'last'),N+2)...
    =-Y(find(xiY>=0.5,1):find(xiY<=1,1,'last'),N+1); % inlet 3
Y(1,:)=Y(2,:); % Left
Y(1,find(yjY>=0.5,1):find(yjY<=1,1,'last'))=...
    2-Y(2,find(yjY>=0.5,1):find(yjY<=1,1,'last')); % inlet 1
Y(M+2,:)=Y(M+1,:); % Right
Y(M+2,find(yjY>=1,1):find(yjY<=1.5,1,'last'))=...
    0.5-Y(M+1,find(yjY>=1,1):find(yjY<=1.5,1,'last')); % inlet 2

%%%%%%%%%% STEP 2 %%%%%%%%%%
[a,b,c,d]=Yvectors(Y,M,N,dt,hx,hy,1/(Re*Sc),2); %Obtain tridiagonal vectors
d=GaussTriSol(a,b,c,d); %Gaussian elimination
% Correspondance with Y, obtain Y(n+1)
for i=2:M+1
    for j=2:N+1
        Y(i,j)=d((i-2)*N+j-1);
    end
end
% Update ghost cells.
Y(:,1)=Y(:,2); % Bottom
Y(:,N+2)=Y(:,N+1); % Top, we will overwrite the inlet 3
Y(find(xiY>=0.5,1):find(xiY<=1,1,'last'),N+2)...
    =-Y(find(xiY>=0.5,1):find(xiY<=1,1,'last'),N+1); % inlet 3
Y(1,:)=Y(2,:); % Left
Y(1,find(yjY>=0.5,1):find(yjY<=1,1,'last'))=...
    2-Y(2,find(yjY>=0.5,1):find(yjY<=1,1,'last')); % inlet 1
Y(M+2,:)=Y(M+1,:); % Right
Y(M+2,find(yjY>=1,1):find(yjY<=1.5,1,'last'))=...
    0.5-Y(M+1,find(yjY>=1,1):find(yjY<=1.5,1,'last')); % inlet 2
end