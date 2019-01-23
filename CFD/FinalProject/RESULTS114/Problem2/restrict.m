function r2h = restrict(rh)
    M=size(rh,1)-2;
    N=size(rh,2)-2;
    M2h=M/2;
    N2h=N/2;
    r2h=zeros(M2h+2,N2h+2);
    rh(:,[1 end])=[];
    rh([1 end],:)=[];
    for i=1:M2h
        for j=1:N2h
            r2h(i+1,j+1) = 0.25*sum(sum(rh(2*i-1:2*i,2*j-1:2*j)));
        end
    end
    % Update the ghost cells
	r2h(1,:) = r2h(2,:);
	r2h(M2h+2,:) = r2h(M2h+1,:);
    r2h(:,1) = r2h(:,2);
	r2h(:,N2h+2) = r2h(:,N2h+1);
end
