function HY=HyperbolicY(Y,M,N,dt,u,v,HY)
Yconv = TVDRK3_2D(Y,M,N,dt,u,v);
HY(4:M+3,4:N+3)=(Yconv(4:M+3,4:N+3)-Y(4:M+3,4:N+3))/dt;
end