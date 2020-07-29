function phi3 = TVDRK3_2D(phi0,M,N,hx,hy,dt,u,v)
% Constants and prealocation
a10=1;
a20=-3/4; a21=1/4;
a30=-1/12; a31=-1/12; a32=2/3;
phi1=zeros(M+6,N+6);
phi2=zeros(M+6,N+6);
phi3=zeros(M+6,N+6);

%%% STEP 1 %%%
for i=4:M+3
    for j=4:N+3
        uij=(u(i-3,j-2)+u(i-2,j-2))/2; % Values of u,v at cell centers
        vij=(v(i-2,j-3)+v(i-2,j-2))/2;
        [DphiDx0,DphiDy0]=WENO5_2D(phi0,uij,vij,i,j,hx,hy); 
        phi1(i,j)=phi0(i,j)-a10*dt*(uij*DphiDx0+vij*DphiDy0);
    end
end
% Update ghost cells
phi1=applyBCs(phi1,M,N,hx,hy,'Y');
%%% STEP 2 %%%
for i=4:M+3
    for j=4:N+3
        uij=(u(i-3,j-2)+u(i-2,j-2))/2; % Values of u,v at cell centers
        vij=(v(i-2,j-3)+v(i-2,j-2))/2;
        [DphiDx0,DphiDy0]=WENO5_2D(phi0,uij,vij,i,j,hx,hy);
        [DphiDx1,DphiDy1]=WENO5_2D(phi1,uij,vij,i,j,hx,hy);
        phi2(i,j)=phi1(i,j)-a20*dt*(uij*DphiDx0+vij*DphiDy0)...
            -a21*dt*(uij*DphiDx1+vij*DphiDy1);
    end
end
% Update ghost cells
phi2=applyBCs(phi2,M,N,hx,hy,'Y');

%%% STEP 3 %%%
for i=4:M+3
    for j=4:N+3
        uij=(u(i-3,j-2)+u(i-2,j-2))/2; % Values of u,v at cell centers
        vij=(v(i-2,j-3)+v(i-2,j-2))/2;
        [DphiDx0,DphiDy0]=WENO5_2D(phi0,uij,vij,i,j,hx,hy);
        [DphiDx1,DphiDy1]=WENO5_2D(phi1,uij,vij,i,j,hx,hy);
        [DphiDx2,DphiDy2]=WENO5_2D(phi2,uij,vij,i,j,hx,hy);
        phi3(i,j)=phi2(i,j)-a30*dt*(uij*DphiDx0+vij*DphiDy0)...
            -a31*dt*(uij*DphiDx1+vij*DphiDy1)...
            -a32*dt*(uij*DphiDx2+vij*DphiDy2);
    end
end
% Update ghost cells
phi3=applyBCs(phi3,M,N,hx,hy,'Y');
end