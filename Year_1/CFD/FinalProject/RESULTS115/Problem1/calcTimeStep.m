function [dt,n]=calcTimeStep(time,dt,tcompare,n,u,v,CFL)
if ( time < tcompare && time+dt >= tcompare )
    dt=tcompare-time;
    n=n+1;
else
    dt=stableTimeStep(u,v,CFL);
end

end