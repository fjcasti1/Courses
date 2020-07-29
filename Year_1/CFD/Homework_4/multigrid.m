function phi = multigrid(phi,f,h)
    M=length(phi)-2;
    phi = GaussSeidel(phi,f,h,1);
    if M>2
        rh = residual(phi,f,h);
        r2h = restrict(rh);
        e2h = zeros(M/2+2,M/2+2);
        e2h = multigrid(e2h,r2h,2*h);
        eh = prolong(e2h);
        phi = phi+eh;
        phi = GaussSeidel(phi,f,h,1);
    end
end