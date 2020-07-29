function phi3 = TVDRK3(phi0,M,hx,dt,t,a)
% Constants and prealocation
a10=1;
a20=-3/4; a21=1/4;
a30=-1/12; a31=-1/12; a32=2/3;
phi1=zeros(M+6,1);
phi2=zeros(M+6,1);
phi3=zeros(M+6,1);

%%% STEP 1 %%%
for i=4:M+3
    phi1(i)=phi0(i)-a10*a(i)*dt*WENO5(phi0,i,hx,a(i));
end
% Update ghost cells
phi1=updateGhostCells(phi1,t,M);

%%% STEP 2 %%%
for i=4:M+3
    phi2(i)=phi1(i)-a20*a(i)*dt*WENO5(phi0,i,hx,a(i))...
        -a21*a(i)*dt*WENO5(phi1,i,hx,a(i));
end
% Update ghost cells
phi2=updateGhostCells(phi2,t,M);

%%% STEP 3 %%%
for i=4:M+3
    phi3(i)=phi2(i)-a30*a(i)*dt*WENO5(phi0,i,hx,a(i))...
        -a31*a(i)*dt*WENO5(phi1,i,hx,a(i))...
        -a32*a(i)*dt*WENO5(phi2,i,hx,a(i));
end
% Update ghost cells
phi3=updateGhostCells(phi3,t+dt,M);
end