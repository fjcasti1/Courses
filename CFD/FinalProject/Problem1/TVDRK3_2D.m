function phi3 = TVDRK3_2D(phi0,M,N,dt,u,v)

% Constants and prealocation
a10=1;
a20=-3/4; a21=1/4;
a30=-1/12; a31=-1/12; a32=2/3;
phi1=zeros(M+6,N+6);
phi2=zeros(M+6,N+6);
phi3=zeros(M+6,N+6);
%%% STEP 1 %%%
[DphiDx0,DphiDy0,uCell,vCell]=WENO5_2D(phi0,u,v,M,N);
for i=4:M+3
    for j=4:N+3
        phi1(i,j)=phi0(i,j)-a10*dt*(uCell(i,j)*DphiDx0(i,j)+vCell(i,j)*DphiDy0(i,j));
    end
end
% Update ghost cells
phi1=applyBCs(phi1,M,N,'Y');
%%% STEP 2 %%%
[DphiDx1,DphiDy1]=WENO5_2D(phi1,u,v,M,N);
for i=4:M+3
    for j=4:N+3
        phi2(i,j)=phi1(i,j)-a20*dt*(uCell(i,j)*DphiDx0(i,j)+vCell(i,j)*DphiDy0(i,j))...
            -a21*dt*(uCell(i,j)*DphiDx1(i,j)+vCell(i,j)*DphiDy1(i,j));
    end
end
% Update ghost cells
phi2=applyBCs(phi2,M,N,'Y');

%%% STEP 3 %%%
[DphiDx2,DphiDy2]=WENO5_2D(phi2,u,v,M,N);
for i=4:M+3
    for j=4:N+3
     phi3(i,j)=phi2(i,j)-a30*dt*(uCell(i,j)*DphiDx0(i,j)+vCell(i,j)*DphiDy0(i,j))...
            -a31*dt*(uCell(i,j)*DphiDx1(i,j)+vCell(i,j)*DphiDy1(i,j))...
            -a32*dt*(uCell(i,j)*DphiDx2(i,j)+vCell(i,j)*DphiDy2(i,j));
    end
end
% Update ghost cells
phi3=applyBCs(phi3,M,N,'Y');
end