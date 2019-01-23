function eh = prolong(e2h)
    M2h=size(e2h,1)-2;
    N2h=size(e2h,2)-2;
    Mh=2*M2h;
    Nh=2*N2h;
    eh=zeros(Mh+2,Nh+2);
    e2h(:,[1 end])=[];
    e2h([1 end],:)=[];
    % Interior
    for i=1:M2h
        for j=1:N2h
            eh(2*i-1+1:2*i+1,2*j-1+1:2*j+1) = e2h(i,j);
        end
    end
    % Calculate the ghost cells
	eh(1,:) = eh(2,:);
	eh(Mh+2,:) = eh(Mh+1,:);
    eh(:,1) = eh(:,2);
	eh(:,Nh+2) = eh(:,Nh+1);
end