function dt=stableTimeStep(u,v,CFL)
global hx;
global hy;
eps=1e-12;
dtu=CFL*min(hx,hy)/(eps+max(max(abs(2*u)))+max(max(abs(v))));
dtv=CFL*min(hx,hy)/(eps+max(max(abs(u)))+max(max(abs(2*v))));
dt=min(dtu,dtv);
end