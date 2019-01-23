function HY=HyperbolicY(Y,M,N,hx,hy,dt,u,v)
Yconv = TVDRK3_2D(Y,M,N,hx,hy,dt,u,v);
HY=zeros(M+6,N+6);
HY(4:M+3,4:N+3)=(Yconv(4:M+3,4:N+3)-Y(4:M+3,4:N+3))/dt;
end