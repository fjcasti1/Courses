function [F,J,y0]=FandJ(N)
x = linspace(0,pi,N+1)';
h=pi/N;
xi = x(2:N);
y0 = xi/pi;

D1 = (gallery('tridiag',N+1,-1,0,1))/(2*h);
D1(1,:) = [];
D1(end,:) = [];
D1(:,1) = [];
D1(:,end) = [];

D2 = (gallery('tridiag',N+1,1,-2,1))/(h^2);
D2(1,:) = [];
D2(end,:) = [];
D2(:,1) = [];
D2(:,end) = [];

% g1 = [zeros(N-2,1)];
% g1(N-2) = (1/(2*h));
% g2 = [zeros(N-2,1)];
% g2(N-2) = (1/(h^2));
% e1 = g2 - g1;

M = D2 - D1;
e1 = [zeros(N-1,1)];
e1(N-1) = ((2-h)/(2*h^2));

F = @(y) (M*y)-cos(y)+e1;
J = @(y) M + diag(sin(y));
end