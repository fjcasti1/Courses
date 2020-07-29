function sigma = sigmaG(alpha,dt,dx)
sigma = 0.5*(EntropyFix(alpha)-(dt/dx)*alpha.^2);
end