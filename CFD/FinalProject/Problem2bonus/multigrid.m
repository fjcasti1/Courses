function phi = multigrid(phi,rhs,hx,hy)
    M=size(phi,1)-2;
    N=size(phi,2)-2;
    phi = GaussSeidel(phi,rhs,hx,hy,1);
    if (M>2 || N>2)
        rh = residual(phi,rhs,hx,hy);
        r2h = restrict(rh);
        e2h = zeros(M/2+2,N/2+2);
        e2h = multigrid(e2h,r2h,2*hx,2*hy);
        eh = prolong(e2h);
        phi = phi+eh;
        phi = GaussSeidel(phi,rhs,hx,hy,1);
    end
end