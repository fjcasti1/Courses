function [phi, LinfR] = SOR (phi,f,h,nIterMax,eps)
    M = size(phi,1)-2;
    N = size(phi,2)-2;
    rho=0.5*(cos(pi/M)+cos(pi/N));
    w=2/(1+sqrt(1-rho^2));
    r = residual(phi,f,h); % residual of the initial guess
    LinfR=zeros(nIterMax,1);
    LinfR(1)=InfNorm(r);
    for n =1:nIterMax
        % Update the interior cells
        for i=2:M+1
            for j=2:N+1
                phi(i,j)=phi(i,j)+w*(0.25*(phi(i+1,j)+phi(i-1,j)...
                    +phi(i,j+1)+phi(i,j-1)-h^2*f(i,j))-phi(i,j));
            end
        end
        % Update the ghost cells
        phi(1,:) = phi(2,:);
        phi(M+2,:) = phi(M+1,:);
        phi(:,1) = phi(:,2);
        phi(:,M+2) = phi(:,M+1);
        % Calculate residuals
        r = residual(phi,f,h);
        LinfR(n+1,1) = InfNorm(r);
        if LinfR(n+1,1)<eps
            break
        end
    end
end