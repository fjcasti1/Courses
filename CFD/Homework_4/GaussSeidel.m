function phi = GaussSeidel (phi,f,h,n)
    for k=1:n
        M=length(phi)-2;
        % Update the interior cells
        for i=2:M+1 
            for j=2:M+1
                phi(i,j)=0.25*(phi(i+1,j)+phi(i-1,j)...
                    +phi(i,j+1)+phi(i,j-1)-h^2*f(i,j));
            end
        end
        % Update the ghost cells
        phi(1,:) = phi(2,:);
        phi(M+2,:) = phi(M+1,:);
        phi(:,1) = phi(:,2);
        phi(:,M+2) = phi(:,M+1);
    end
end