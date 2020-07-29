function [phi,LinfR] = poisson (phi,rhs)
global hx; global hy;
global nIterMax;
global tol;

r = residual(phi,rhs,hx,hy);
LinfR=zeros(nIterMax,1);
LinfR(1)=InfNorm(r);
for n=1:nIterMax
    phi = multigrid(phi,rhs,hx,hy);
    r = residual(phi,rhs,hx,hy);
    LinfR(n+1,1) = InfNorm(r);
    if LinfR(n+1,1)<tol
        break
    end
end
LinfR(n+2:nIterMax)=[];
end