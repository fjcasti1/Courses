function G=HY_G(u,alpha,dt,dx,M)
sigma=sigmaG(alpha,dt,dx);
G=zeros(M+4,1);
for i=3:M+2
    S=sign(u(i+1)-u(i));
    G(i) = S*max(0,min(sigma(i)*abs(u(i+1)-u(i)),...
            S*sigma(i-1)*(u(i)-u(i-1))));
end
% Periodic Boundary Conditions
G(1)=G(M+1);
G(2)=G(M+2);
G(M+3)=G(3);
G(M+4)=G(4);
end