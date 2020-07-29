function alpha = HY_alpha(u,E,M)
eps=0.1;
alpha = zeros(M+3,1);
for i=1:M+3
    absDelta=abs(u(i+1)-u(i));
    if absDelta>=eps
        alpha(i) = (E(i+1)-E(i))/(u(i+1)-u(i));
    elseif absDelta<eps
        alpha(i) =(u(i)+u(i+1))/2;
    end
end