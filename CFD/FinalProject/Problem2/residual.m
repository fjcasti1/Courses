function r = residual(phi,rhs,hx,hy)
    M = size(phi,1)-2;
    N = size(phi,2)-2;
    r= zeros(size(phi));
    % Calculate the interior cells
    for i=2:M+1
        for j=2:N+1
            r(i,j)=rhs(i,j)-(phi(i+1,j)-2*phi(i,j)+phi(i-1,j))/hx^2 ...
            -(phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/hy^2;
        end
    end
    % calculate the ghost cells
    r(1,:) = r(2,:);
    r(M+2,:) = r(M+1,:);
    r(:,1) = r(:,2);
    r(:,N+2) = r(:,N+1);
end