function [phi,LinfR] = poisson (phi,f,h,nIterMax,eps)
    [M,N] = size(phi);
    r = residual(phi,f,h);
    LinfR=zeros(nIterMax,1);
    LinfR(1)=InfNorm(r);
    for n=1:nIterMax
        phi = multigrid(phi,f,h);
        r = residual(phi,f,h);
        LinfR(n+1,1) = InfNorm(r);
        if LinfR(n+1,1)<eps
            break
        end
    end
    LinfR(n+2:nIterMax)=[];
end