function u = TVD2order(u,h,dt,dx,M)
for i=3:M+2
   u(i)=u(i)- dt/dx*(h(i-1)-h(i-2));
end
% Periodic Boundary Conditions
u(1)=u(M+1);
u(2)=u(M+2);
u(M+3)=u(3);
u(M+4)=u(4);
end