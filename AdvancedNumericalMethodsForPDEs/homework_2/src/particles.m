N=10000;

% compute N particles with p[xi=x]=x/8 for 0<x<2
eta =rand (1,N);
xi=sqrt(4*eta);

figure(1)
plot(xi,rand(1,N),'.')
xlabel('xi')


M=100;
dx=1/(M-1);
%compute u(x) with M-1 cells
x=linspace(0,2,M);


for m=1:M-1
    I=find((x(m)<xi)&(xi<x(m+1)));
    u(m)=length(I)/N/dx;
end
figure(2)
plot(x(1:M-1),u,'r.')
xlabel('x')
ylabel('u')