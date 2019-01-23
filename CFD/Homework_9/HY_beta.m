function beta = HY_beta(u,G,M)
eps=1e-12;
beta=zeros(M+1,1);
for i=2:M+2
    absDelta=abs(u(i+1)-u(i));
    if absDelta>=eps
        beta(i-1) = (G(i+1)-G(i))/(u(i+1)-u(i));
    elseif absDelta<eps
        beta(i-1)=0;
    end
end
end