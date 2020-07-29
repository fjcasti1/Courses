function Phi = HY_Phi(u,G,alpha,beta,M)
psi=EntropyFix(alpha(2:end-1)+beta);
Phi=zeros(M+1,1);
for i=2:M+2
    Phi(i-1)=G(i+1)+G(i)-psi(i-1)*(u(i+1)-u(i));
end
end