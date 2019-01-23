function [u,v,Hu,Hv]=solveBurgers2D(u,v,Hu,Hv,M,N,dt,time)
Huold=Hu;
Hvold=Hv;
[Hu,Hv]=HyperbolicBurgers2D(u,v,M,N);
if time==0 %FCTS for the first time step
    Huold=Hu;
    Hvold=Hv;
end
u=ADI_u(u,M,N,dt,time,1.5*Hu-0.5*Huold);
v=ADI_v(v,M,N,dt,time,1.5*Hv-0.5*Hvold);

end