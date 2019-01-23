function r = residual(phi,f,h)
    M = size(phi,1)-2;
%     M
    r= zeros(size(phi));
%     whos phi f r
%     keyboard
    % Calculate the interior cells
    for i=2:M+1
        for j=2:M+1
            r(i,j)=f(i,j)-(phi(i+1,j)-2*phi(i,j)+phi(i-1,j)+phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/h^2;
        end
    end
    % calculate the ghost cells
    r(1,:) = r(2,:);
    r(M+2,:) = r(M+1,:);
    r(:,1) = r(:,2);
    r(:,M+2) = r(:,M+1);
end