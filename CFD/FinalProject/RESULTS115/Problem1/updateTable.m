function [T]=updateTable(T)
% T is a matrix that will be turned into a table (outside) this function
phi=T(:,3);
r=2;
Fsec=1.25;

% GCI
for i=3:length(phi)
    p=log(abs(phi(i-2)-phi(i-1))/abs(phi(i-1)-phi(i)))/log(r);
    phi_h0=phi(i)+(phi(i)-phi(i-1))/(r^p-1);
    GCI12=Fsec*abs(1-phi(i-1)/phi(i))/(r^p-1)*100;
    GCI23=Fsec*abs(1-phi(i-2)/phi(i-1))/(r^p-1)*100;
    beta=GCI12*r^p/GCI23;
    % Include results in a table
    T(i,4)=p;
    T(i,5)=phi_h0;
    T(i,6)=GCI12;
    T(i,7)=GCI23;
    T(i,8)=beta;
end

end