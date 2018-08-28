function basin = basin_henon(n,x,y,param)
basin  = zeros(n,n);
[X,Y]  = meshgrid(x,y);
v(:,1) = reshape(X,[n^2,1]);
v(:,2) = reshape(Y,[n^2,1]);
% param.maxiter=11
for k=1:param.maxiter
    v = henon(v,param);
     if(ismember(1,isnan(v)))
%          keyboard
     end 
end
% v
% keyboard
basin = double(reshape(abs((v(:,1))>param.L | abs(v(:,2))>param.L),[n,n]));
% keyboard
end