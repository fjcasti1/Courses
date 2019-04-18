close all;
N=40;
[D,x] = cheb(N);
Dbc=D;
Dbc([1 end],:)=0;
NBC = zeros(size(x));
NBC([1 end])=1;
M=D*Dbc+diag(x.^2)*Dbc+eye(size(x));
f=-(D*NBC-diag(x.^2)*NBC);
u=M\f;
figure
plot(x,u)