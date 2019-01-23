function phi = GaussSeidel (phi,rhs,hx,hy,n)
    for k=1:n
        M=size(phi,1)-2;
        N=size(phi,2)-2;
        h1=hy^2/(2*(hx^2+hy^2));
        h2=hx^2/(2*(hx^2+hy^2));
        h3=hx^2*hy^2/(2*(hx^2+hy^2));
        % Update the interior cells
        for i=2:M+1 
            for j=2:N+1
                phi(i,j)=h1*(phi(i+1,j)+phi(i-1,j))...
                    +h2*(phi(i,j+1)+phi(i,j-1))-h3*rhs(i,j);
            end
        end
        % Update the ghost cells
        phi(1,:) = phi(2,:);
        phi(M+2,:) = phi(M+1,:);
        phi(:,1) = phi(:,2);
        phi(:,N+2) = phi(:,N+1);
    end
end